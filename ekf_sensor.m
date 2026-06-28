% ekf_sensor.m
% Step 5a — Extended Kalman Filter (EKF) soft sensor
% Fuses Raman PLS predictions with Monod kinetics for optimal state estimation
% Adapted to indpensim_run() data structure (batch_data{i}.X.y, .S.y, etc.)

project_dir = 'C:\Users\suraj\Documents\NITS\6th Sem\SNBose\bioreactor_project';
cd(project_dir);
addpath('C:\Users\suraj\Documents\MATLAB\IndPenSim_V2.01');
load(fullfile(project_dir, 'batch_dataset.mat'));   % batch_data, N
load(fullfile(project_dir, 'soft_sensor.mat'));     % beta_final, nComp, nWave

set(0, 'DefaultFigureWindowStyle', 'normal');

%% --- Parameters ---
% Use batch 81 (first test batch) for EKF demonstration
test_batch_idx = 81;
batch = batch_data{test_batch_idx};

% Time vector
t_vec = batch.X.t;           % [1150 x 1]
dt    = mean(diff(t_vec));    % average time step (h)
T_meas = numel(t_vec);

fprintf('Running EKF on batch %d  |  %d time steps  |  dt = %.3f h\n', ...
    test_batch_idx, T_meas, dt);

%% --- Monod kinetic parameters (nominal IndPenSim values) ---
mu_max = 0.092;   % max specific growth rate (1/h)
Ks     = 0.100;   % Monod saturation constant (g/L)
Yxs    = 0.450;   % biomass yield on substrate (g/g)
qp_max = 0.005;   % max specific penicillin production rate (1/h)
Kp     = 0.020;   % penicillin saturation constant (g/L)

%% --- EKF noise tuning ---
% Process noise Q: how much we trust the kinetic model
% Measurement noise R: how much we trust the Raman PLS sensor
Q       = diag([0.50, 0.50, 1e-4]);   % high process noise = trust sensor more
R_noise = 1.85^2;                      % Raman PLS variance = RMSEP^2 from Step 4

%% --- Initialise state and covariance ---
% State vector: z = [X (biomass), S (substrate), P (penicillin)]
z_hat = [batch.X.y(1);   % initial biomass
         batch.S.y(1);    % initial substrate
         batch.P.y(1)];   % initial penicillin
P_cov = eye(3) * 0.1;

X_ekf = zeros(T_meas, 1);   % EKF biomass estimates
S_ekf = zeros(T_meas, 1);   % EKF substrate estimates
P_ekf = zeros(T_meas, 1);   % EKF penicillin estimates

%% --- Prepare Raman measurements for this batch ---
% Transpose Intensity to [nTime x nWave], then SNV-normalise
R_raw = batch.Raman_Spec.Intensity';   % [1150 x 2200]
R_mean_row = mean(R_raw, 2);
R_std_row  = std(R_raw, 0, 2);
zero_rows  = R_std_row < 1e-10;
R_std_row(zero_rows) = 1;
R_snv = (R_raw - R_mean_row) ./ R_std_row;
R_snv(zero_rows, :) = 0;

% PLS prediction: Raman -> biomass
X_raman = [ones(T_meas, 1), R_snv] * beta_final;   % [1150 x 1]

fprintf('Raman PLS prediction range: [%.2f, %.2f] g/L\n', ...
    min(X_raman), max(X_raman));

%% --- EKF Loop ---
for k = 1:T_meas

    % --- Inputs at this time step ---
    Fs_k = batch.Fs.y(k);   % substrate feed rate (L/h)
    V_k  = batch.V.y(k);    % volume (L) — needed for dilution term

    % Guard against zero volume
    if V_k < 1e-6; V_k = 1e-6; end

    % --- Monod kinetics ---
    mu   = mu_max * z_hat(2) / (Ks  + z_hat(2));
    q_p  = qp_max * z_hat(1) / (Kp  + z_hat(1));

    % --- Process model (Euler integration) ---
    dX = mu  * z_hat(1);
    dS = -dX / Yxs + Fs_k / V_k;
    dP = q_p;

    z_pred = max(z_hat + dt * [dX; dS; dP], 0);   % enforce non-negative

    % --- Jacobian F = dF/dz (linearised process model) ---
    dmu_dS  = mu_max * Ks / (Ks + z_hat(2))^2;
    dqp_dX  = qp_max * Kp / (Kp + z_hat(1))^2;

    F_jac = eye(3) + dt * [
        mu,                  dmu_dS * z_hat(1),   0;
       -mu / Yxs,           -dmu_dS * z_hat(1) / Yxs,  0;
        dqp_dX,              0,                   0 ];

    % --- Predict covariance ---
    P_pred = F_jac * P_cov * F_jac' + Q;

    % --- Measurement update: observe biomass from Raman PLS ---
    H     = [1, 0, 0];                          % observe X only
    z_meas = X_raman(k);                        % Raman PLS estimate
    S_innov = H * P_pred * H' + R_noise;        % innovation covariance
    K_g   = P_pred * H' / S_innov;             % Kalman gain

    z_hat = z_pred + K_g * (z_meas - H * z_pred);
    z_hat = max(z_hat, 0);                      % enforce non-negative
    P_cov = (eye(3) - K_g * H) * P_pred;

    % --- Store estimates ---
    X_ekf(k) = z_hat(1);
    S_ekf(k) = z_hat(2);
    P_ekf(k) = z_hat(3);
end

%% --- Evaluate EKF vs raw Raman PLS ---
X_actual = batch.X.y(1:T_meas);

RMSE_raman = sqrt(mean((X_raman   - X_actual).^2));
RMSE_ekf   = sqrt(mean((X_ekf     - X_actual).^2));
R2_raman   = 1 - sum((X_raman - X_actual).^2) / sum((X_actual - mean(X_actual)).^2);
R2_ekf     = 1 - sum((X_ekf   - X_actual).^2) / sum((X_actual - mean(X_actual)).^2);

fprintf('\n--- EKF vs Raw Raman PLS (Batch %d) ---\n', test_batch_idx);
fprintf('%-20s  RMSE = %.4f g/L   R2 = %.4f\n', 'Raman PLS only:', RMSE_raman, R2_raman);
fprintf('%-20s  RMSE = %.4f g/L   R2 = %.4f\n', 'EKF fused:',     RMSE_ekf,   R2_ekf);
if RMSE_ekf < RMSE_raman
    fprintf('EKF improved RMSE by %.1f%%\n', (RMSE_raman - RMSE_ekf)/RMSE_raman*100);
else
    fprintf('Note: EKF did not improve RMSE — consider tuning Q/R_noise\n');
end

%% --- Plot 1: Biomass — Actual vs Raman PLS vs EKF ---
figure;
plot(t_vec, X_actual, 'b-',  'LineWidth', 2.0); hold on;
plot(t_vec, X_raman,  'r--', 'LineWidth', 1.5);
plot(t_vec, X_ekf,    'g-',  'LineWidth', 1.5);
xlabel('Time (h)');
ylabel('Biomass X (g/L)');
title(sprintf('Batch %d — Biomass: Actual vs Raman PLS vs EKF', test_batch_idx));
legend('Actual', 'Raman PLS', 'EKF Estimate', 'Location', 'southeast');
grid on;

%% --- Plot 2: EKF substrate and penicillin estimates ---
figure;
subplot(2,1,1);
plot(t_vec, batch.S.y(1:T_meas), 'b-', 'LineWidth', 1.5); hold on;
plot(t_vec, S_ekf, 'r--', 'LineWidth', 1.5);
ylabel('Substrate S (g/L)');
title(sprintf('Batch %d — EKF State Estimates', test_batch_idx));
legend('Actual S', 'EKF S');
grid on;

subplot(2,1,2);
plot(t_vec, batch.P.y(1:T_meas), 'b-', 'LineWidth', 1.5); hold on;
plot(t_vec, P_ekf, 'r--', 'LineWidth', 1.5);
xlabel('Time (h)');
ylabel('Penicillin P (g/L)');
legend('Actual P', 'EKF P');
grid on;

%% --- Save EKF results ---
save(fullfile(project_dir, 'ekf_results.mat'), ...
    'X_ekf', 'S_ekf', 'P_ekf', 'X_raman', 'X_actual', ...
    'RMSE_ekf', 'RMSE_raman', 'R2_ekf', 'R2_raman', ...
    'test_batch_idx', '-v7.3');
fprintf('\nSaved to ekf_results.mat\n');