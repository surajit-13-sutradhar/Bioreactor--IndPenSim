% soft_sensor_pls.m
% PLS soft sensor for biomass prediction from Raman spectra
% Raman_Spec is a struct with:
%   .Wavelength  [2200 x 1]    — wavenumber axis
%   .Intensity   [2200 x 1150] — rows = time points, cols = wavenumbers

project_dir = 'C:\Users\suraj\Documents\NITS\6th Sem\SNBose\bioreactor_project';
cd(project_dir);
addpath('C:\Users\suraj\Documents\MATLAB\IndPenSim_V2.01');
load(fullfile(project_dir, 'batch_dataset.mat'));  % loads batch_data, N

set(0, 'DefaultFigureWindowStyle', 'normal');

%% --- Step 1: Confirm Raman structure ---
R_sample   = batch_data{1}.Raman_Spec.Intensity';  % transpose → [1150 x 2200]
nTime      = size(R_sample, 1);                     % 1150 time points
nWave      = size(R_sample, 2);                     % 2200 wavenumbers
wavelength = batch_data{1}.Raman_Spec.Wavelength;   % [2200 x 1]

fprintf('Raman Intensity size after transpose: [%d x %d]  (time x wavenumbers)\n', nTime, nWave);

% Sanity check: nTime must match biomass time points
nTime_X = numel(batch_data{1}.X.y);
if nTime ~= nTime_X
    warning('Raman nTime (%d) != Biomass nTime (%d) — will use min.', nTime, nTime_X);
    nTime = min(nTime, nTime_X);
end

%% --- Step 2: Stack all spectra and biomass values ---
R_all = zeros(N * nTime, nWave);
X_all = zeros(N * nTime, 1);

for i = 1:N
    rows        = (i-1)*nTime + (1:nTime);
    R_all(rows, :) = batch_data{i}.Raman_Spec.Intensity(1:nWave, 1:nTime)'; % transpose
    X_all(rows)    = batch_data{i}.X.y(1:nTime);
end

fprintf('Stacked spectra:  [%d x %d]\n', size(R_all));
fprintf('Stacked biomass:  [%d x 1]\n',  size(X_all, 1));

%% --- Step 3: SNV normalisation (row-wise: per spectrum) ---
R_mean = mean(R_all, 2);
R_std  = std(R_all,  0, 2);

% Replace zero-std rows (flat spectra) with zeros to avoid NaN
zero_rows       = R_std < 1e-10;
R_std(zero_rows) = 1;
R_snv           = (R_all - R_mean) ./ R_std;
R_snv(zero_rows, :) = 0;

if any(zero_rows)
    fprintf('Warning: %d flat-spectrum rows set to zero after SNV.\n', sum(zero_rows));
end

%% --- Step 4: Train/test split (80 batches train, 20 test) ---
n_train = round(0.8 * N);
train_end = n_train * nTime;

R_train = R_snv(1:train_end, :);      Y_train = X_all(1:train_end);
R_test  = R_snv(train_end+1:end, :);  Y_test  = X_all(train_end+1:end);

fprintf('\nTrain: %d samples (%d batches)\n', size(R_train,1), n_train);
fprintf('Test:  %d samples (%d batches)\n',  size(R_test,1),  N-n_train);

%% --- Step 5: Cross-validate number of PLS components ---
fprintf('\nRunning 10-fold CV for component selection (this may take a minute)...\n');
max_comp = 12;
RMSECV   = zeros(max_comp, 1);

for nc = 1:max_comp
    [~,~,~,~,~, MSECV] = plsregress(R_train, Y_train, nc, 'CV', 10);
    RMSECV(nc) = sqrt(MSECV(end));  % MSECV(end) = prediction error at nc components
end

[~, best_nComp] = min(RMSECV);
fprintf('Best PLS components (CV): %d\n', best_nComp);

figure;
plot(1:max_comp, RMSECV, 'bo-', 'LineWidth', 1.5);
xlabel('Number of PLS Components');
ylabel('RMSECV (g/L)');
title('Cross-Validation: PLS Component Selection');
xline(best_nComp, 'r--', sprintf('Best = %d', best_nComp));
grid on;

%% --- Step 6: Fit final PLS model ---
nComp = best_nComp;
[~,~,~,~, beta_final] = plsregress(R_train, Y_train, nComp);

%% --- Step 7: Predict and evaluate ---
Y_pred_train = [ones(size(R_train,1),1), R_train] * beta_final;
Y_pred_test  = [ones(size(R_test,1),1),  R_test]  * beta_final;

RMSEC = sqrt(mean((Y_pred_train - Y_train).^2));
RMSEP = sqrt(mean((Y_pred_test  - Y_test).^2));
R2_tr = 1 - sum((Y_pred_train - Y_train).^2) / sum((Y_train - mean(Y_train)).^2);
R2    = 1 - sum((Y_pred_test  - Y_test).^2)  / sum((Y_test  - mean(Y_test)).^2);

fprintf('\n--- PLS Soft Sensor Results ---\n');
fprintf('Training:  RMSEC = %.4f g/L,  R2 = %.4f\n', RMSEC, R2_tr);
fprintf('Test:      RMSEP = %.4f g/L,  R2 = %.4f\n', RMSEP, R2);
if R2 > 0.90
    fprintf('Model quality: GOOD (R2 > 0.90)\n');
elseif R2 > 0.80
    fprintf('Model quality: ACCEPTABLE (R2 > 0.80)\n');
else
    fprintf('Model quality: POOR — consider more components or check spectra\n');
end

%% --- Step 8: Scatter plot ---
figure;
scatter(Y_test, Y_pred_test, 20, 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
ref_range = [min(Y_test), max(Y_test)];
plot(ref_range, ref_range, 'r-', 'LineWidth', 1.5);
xlabel('Actual Biomass (g/L)');
ylabel('Predicted Biomass (g/L)');
title(sprintf('PLS Soft Sensor — R^2 = %.3f,  RMSEP = %.3f g/L', R2, RMSEP));
legend('Test predictions', 'Perfect fit', 'Location', 'northwest');
grid on;

%% --- Step 9: Time-series plot for first test batch ---
figure;
t_axis = batch_data{1}.X.t(1:nTime);
plot(t_axis, Y_test(1:nTime),      'b-',  'LineWidth', 1.5); hold on;
plot(t_axis, Y_pred_test(1:nTime), 'r--', 'LineWidth', 1.5);
xlabel('Time (h)');
ylabel('Biomass X (g/L)');
title(sprintf('Batch %d — Actual vs PLS-Predicted Biomass', n_train+1));
legend('Actual', 'PLS Predicted');
grid on;

%% --- Save ---
save(fullfile(project_dir, 'soft_sensor.mat'), ...
    'beta_final', 'nComp', 'nWave', 'wavelength', 'RMSEP', 'R2', '-v7.3');
fprintf('\nSaved to soft_sensor.mat\n');