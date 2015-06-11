function effiency_30C()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get and clean the data for both 30x30 and Chimera plots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 30 x 30
load_file = ['experiment_2_n' int2str(100) '.mat'];
grid_data = load(load_file);
grid_data = grid_data.recorder;

for i = [100 1000 10000]
    load_file = ['experiment_2_n' int2str(i) '.mat'];
    current_grid_data = load(load_file);
    grid_data = [grid_data ; current_grid_data.recorder];
end

%% Chimera 
load_file = ['experiment_4_n' int2str(100) '.mat'];
C_data = load(load_file);
C_data = C_data.recorder;

for i = [100 1000 10000]
    load_file = ['experiment_4_n' int2str(i) '.mat'];
    current_C_data = load(load_file);
    C_data = [C_data ; current_C_data.recorder];
end

%%%%%%%%%%%%%%%%%%%%%%%
% Construct the plots %
%%%%%%%%%%%%%%%%%%%%%%%
%% 30 x 30

%% Initialize data-recording variables
%% recorder is (ntrials, 7) matrix:
%% recorder(:,1) = sample size
%% recorder(:,2) = mle factorization size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf factorization size
%% recorder(:,6) = mpf sq error
%% recorder(:,7) = mpf time
sizes = unique(grid_data(:,1));

grid_mean_results = nan(size(sizes, 1), size(grid_data, 2));
grid_std_results = nan(size(sizes, 1), size(grid_data, 2));

for i = 1:length(sizes)
    current_grid_data = grid_data(grid_data(:,1) == sizes(i),:);
    grid_mean_results(i, :) = mean(current_grid_data);
    grid_std_results(i, :) = std(current_grid_data);
end

%% Chimera
%% Initialize data-recording variables
%% recorder is (ntrials, 7) matrix:
%% recorder(:,1) = sample size
%% recorder(:,2) = mle factorization size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf factorization size
%% recorder(:,6) = mpf sq error
%% recorder(:,7) = mpf time
sizes = unique(C_data(:,1));

C_mean_results = nan(size(sizes, 1), size(C_data, 2));
C_std_results = nan(size(sizes, 1), size(C_data, 2));

for i = 1:length(sizes)
    current_C_data = C_data(C_data(:,1) == sizes(i),:);
    C_mean_results(i, :) = mean(current_C_data);
    C_std_results(i, :) = std(current_C_data);
end

clf;
hold on;
grid on;
set(gca, 'xscale', 'log', 'yscale', 'log');
xlim([10^1.9, 10^4.1]);
ylim([10^-0.7, 10^2.2]);
axis square;

%% 4 x 4
errorbar(grid_mean_results(:, 1), grid_mean_results(:, 3), grid_std_results(:, 3), ...
         'o-', 'color', [0 0 0], 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         [0 0 0]);
errorbar(grid_mean_results(:, 1), grid_mean_results(:, 6), grid_std_results(:, 6), ...
         'o-', 'color', [0.5 0.5 0.2], 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         [0.5 0.5 0.2]);
%% Complete graph
errorbar(C_mean_results(:, 1), C_mean_results(:, 3), C_std_results(:, 3), ...
         'ro-', 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         'r');
errorbar(C_mean_results(:, 1), C_mean_results(:, 6), C_std_results(:, 6), ...
         'bo-', 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         'b');
xlabel('Sample size', 'interpreter', 'latex', 'fontsize', 25)
ylabel('Mean Squared Error', 'interpreter', 'latex', 'fontsize', 25)
h = legend('MLE Grid', 'MPF Grid', ...
           'MLE Chimera', 'MPF Chimera', 'location', 'northeast');
set(h, 'interpreter', 'latex', 'fontsize', 25);


%% Save figure
saveas(gcf,'plots/effiency_30C','epsc')