% run_baseline.m
addpath('C:\Users\suraj\Documents\MATLAB\IndPenSim_V2.01');

% Build the Batch_run_flags struct that indpensim_run expects
Batch_run_flags.Control_strategy          = 0;  % 0 = recipe driven (no PRBS)
Batch_run_flags.Batch_length              = 0;  % 0 = fixed batch length (230h)
Batch_run_flags.Batch_fault_order_reference = 0;  % 0 = no faults
Batch_run_flags.Raman_spec                = 1;  % 1 = record Raman spectra

% Run batch number 1
Batch_no = 1;
batch = indpensim_run(Batch_no, Batch_run_flags);

% Override the docked figure setting from indpensim_run
set(0, 'DefaultFigureWindowStyle', 'normal');

% Plot penicillin concentration over time
figure;
plot(batch.P.t, batch.P.y, 'b', 'LineWidth', 1.5);
xlabel('Time (h)');
ylabel('Penicillin (g/L)');
title('Baseline Batch - Penicillin Profile');
grid on;

% Print yield statistics
disp(batch.Stats);

% Show all available field names for use in later steps
disp(fieldnames(batch));