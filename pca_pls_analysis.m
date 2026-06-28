% pca_pls_analysis.m
project_dir = 'C:\Users\suraj\Documents\NITS\6th Sem\SNBose\bioreactor_project';
cd(project_dir);
addpath('C:\Users\suraj\Documents\MATLAB\IndPenSim_V2.01');
load(fullfile(project_dir, 'eda_results.mat'));

set(0, 'DefaultFigureWindowStyle', 'normal');

% Remove RPM (constant, NaN correlation) from analysis
vars_clean = {'Fs','Fg','Fb','Fc','S','DO2','X','pH','T','PAA','NH3','CO2outgas','OUR','CER'};
rpm_idx    = strcmp(vars, 'RPM');
X3_clean   = X3(:,:, ~rpm_idx);       % remove RPM slice
X_mean_clean = X_mean(:, ~rpm_idx);   % remove RPM column
nVars_clean  = numel(vars_clean);

%% ---- PCA ----
% Unfold 3D array to 2D and standardise
X_unfold = reshape(X3_clean, N, nTime * nVars_clean);
X_std    = zscore(X_unfold);

[coeff, score, ~, ~, explained] = pca(X_std);

% Scree plot
figure;
subplot(1,2,1);
bar(cumsum(explained(1:15)));
xlabel('Principal Component');
ylabel('Cumulative Variance Explained (%)');
title('PCA Scree Plot');
yline(90, 'r--', '90%');
grid on;

% Score plot coloured by yield
subplot(1,2,2);
scatter(score(:,1), score(:,2), 60, P_final, 'filled');
colorbar;
xlabel('PC1'); ylabel('PC2');
title('PCA Score Plot – Coloured by Penicillin Yield');
grid on;

% Print variance explained by first 5 PCs
fprintf('\nPCA – Variance explained:\n');
fprintf('%-5s  %-20s  %-25s\n', 'PC', 'Individual (%)', 'Cumulative (%)');
fprintf('%s\n', repmat('-', 1, 45));
for k = 1:5
    fprintf('%-5d  %-20.2f  %-25.2f\n', k, explained(k), sum(explained(1:k)));
end

%% ---- PLS ----
X_pls = zscore(X_mean_clean);
Y_pls = zscore(P_final);

[~,~,~,~, BETA, PCTVAR,~, stats] = plsregress(X_pls, Y_pls, 5);

% VIP scores
W0    = stats.W ./ sqrt(sum(stats.W.^2));
p     = size(X_pls, 2);
sumSq = sum(diag(X_pls' * X_pls) .* BETA(2:end,:).^2);
VIP   = sqrt(p * sum(W0.^2 .* sumSq, 2) / sum(sumSq));

% VIP bar chart
figure;
bar(VIP);
set(gca, 'XTickLabel', vars_clean, 'XTickLabelRotation', 45);
yline(1, 'r--', 'VIP = 1');
ylabel('VIP Score');
title('PLS VIP – CPP Importance for Penicillin Yield');
grid on;

% Print VIP table
fprintf('\nPLS VIP Scores (sorted):\n');
fprintf('%-15s  VIP\n', 'Variable');
fprintf('%s\n', repmat('-', 1, 25));
[vip_sorted, vip_idx] = sort(VIP, 'descend');
for k = 1:nVars_clean
    flag = '';
    if vip_sorted(k) > 1; flag = ' <-- confirmed CPP'; end
    fprintf('%-15s  %.4f%s\n', vars_clean{vip_idx(k)}, vip_sorted(k), flag);
end

save(fullfile(project_dir, 'pca_pls_results.mat'), ...
    'coeff', 'score', 'explained', 'VIP', 'vars_clean', 'X_pls', 'Y_pls', '-v7.3');
disp('Saved to pca_pls_results.mat');