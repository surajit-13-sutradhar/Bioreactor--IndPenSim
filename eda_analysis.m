% eda_analysis.m
project_dir = 'C:\Users\suraj\Documents\NITS\6th Sem\SNBose\bioreactor_project';
cd(project_dir);
addpath('C:\Users\suraj\Documents\MATLAB\IndPenSim_V2.01');
load(fullfile(project_dir, 'batch_dataset.mat'));

set(0, 'DefaultFigureWindowStyle', 'normal');

vars  = {'Fs','Fg','RPM','Fb','Fc','S','DO2','X','pH','T','PAA','NH3','CO2outgas','OUR','CER'};
nVars = numel(vars);
nTime = numel(batch_data{1}.P.t);

% Build 3D array [batches x time x variables]
X3 = zeros(N, nTime, nVars);
for i = 1:N
    for v = 1:nVars
        X3(i,:,v) = batch_data{i}.(vars{v}).y;
    end
end

% Final penicillin yield per batch
P_final = arrayfun(@(i) batch_data{i}.P.y(end), 1:N)';

% Time-averaged values then Pearson correlation with yield
X_mean    = squeeze(mean(X3, 2));
corr_vals = corr(X_mean, P_final);

% Ranked bar chart
[sorted_corr, idx] = sort(abs(corr_vals), 'descend');
figure;
bar(abs(corr_vals(idx)));
set(gca, 'XTickLabel', vars(idx), 'XTickLabelRotation', 45);
ylabel('|Pearson r| with Final Penicillin Yield');
title('CPP Ranking by Correlation with CQA (Penicillin Yield)');
yline(0.5, 'r--', 'Threshold = 0.5');
grid on;

% Print ranked table
fprintf('\nCPP Ranking (sorted by correlation with yield):\n');
fprintf('%-15s  |r|\n', 'Variable');
fprintf('%s\n', repmat('-', 1, 25));
for k = 1:nVars
    fprintf('%-15s  %.4f\n', vars{idx(k)}, sorted_corr(k));
end

save(fullfile(project_dir, 'eda_results.mat'), 'X3', 'X_mean', 'P_final', 'vars', 'nVars', 'nTime', 'N', '-v7.3');
disp('Saved to eda_results.mat');