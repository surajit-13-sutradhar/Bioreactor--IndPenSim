% generate_batches.m
addpath('C:\Users\suraj\Documents\MATLAB\IndPenSim_V2.01');

N = 100;
batch_data = cell(N, 1);

% Same flags for all batches - no faults, record Raman spectra
Batch_run_flags.Control_strategy            = zeros(1, N);  % 0 = recipe driven
Batch_run_flags.Batch_length                = zeros(1, N);  % 0 = fixed 230h
Batch_run_flags.Batch_fault_order_reference = zeros(1, N);  % 0 = no faults
Batch_run_flags.Raman_spec                  = ones(1, N);   % 1 = record Raman

for i = 1:N
    fprintf('Running batch %d of %d...\n', i, N);
    batch_data{i} = indpensim_run(i, Batch_run_flags);
end

set(0, 'DefaultFigureWindowStyle', 'normal');

% Plot all penicillin profiles overlaid to check variability
figure;
for i = 1:N
    plot(batch_data{i}.P.t, batch_data{i}.P.y, 'Color', [0.5 0.5 0.5 0.3]);
    hold on;
end
xlabel('Time (h)');
ylabel('Penicillin (g/L)');
title('All 100 Batches – Penicillin Profiles');
grid on;

% Print yield summary across all batches
yields = arrayfun(@(i) batch_data{i}.P.y(end), 1:N);
fprintf('\nYield across 100 batches:\n');
fprintf('  Mean:  %.4f g/L\n', mean(yields));
fprintf('  Std:   %.4f g/L\n', std(yields));
fprintf('  Min:   %.4f g/L\n', min(yields));
fprintf('  Max:   %.4f g/L\n', max(yields));

save('batch_dataset.mat', 'batch_data', 'N');
disp('Saved to batch_dataset.mat');