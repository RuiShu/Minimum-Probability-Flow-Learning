function effiency_4K()
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get and clean the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4 by 4 grid
load_file = ['experiment_1_n' int2str(100) '.mat'];
grid_data = load(load_file);
grid_data = grid_data.recorder;

for i = [1000 10000]
    load_file = ['experiment_1_n' int2str(i) '.mat'];
    current_grid_data = load(load_file);
    grid_data = [grid_data ; current_grid_data.recorder];
end

%% K graphs 
load_file = ['experiment_3_n' int2str(100) '_nodes6.mat'];
K_data = load(load_file);
K_data = K_data.recorder;

for i = [1000 10000]
    load_file = ['experiment_3_n' int2str(i) '_nodes6.mat'];
    current_K_data = load(load_file);
    K_data = [K_data ; current_K_data.recorder];
end

%%%%%%%%%%%%%%%%%%%
% Plot the curves %
%%%%%%%%%%%%%%%%%%%
%% Initialize data-recording variables
%% recorder is (ntrials, 5) matrix:
%% recorder(:,1) = nsamples
%% recorder(:,2) = mle sq error
%% recorder(:,3) = mle time
%% recorder(:,4) = mpf sq error
%% recorder(:,5) = mpf time
sizes = unique(grid_data(:,1));

meanresults = nan(size(sizes, 1), 5);
stdresults = nan(size(sizes, 1), 5);

for i = 1:length(sizes)
    current_grid_data = grid_data(grid_data(:,1) == sizes(i),:);
    grid_mean_results(i, :) = mean(current_grid_data);
    grid_std_results(i, :) = std(current_grid_data);
end

%% Initialize data-recording variables
%% recorder is (ntrials, 6) matrix:
%% recorder(:,1) = nsamples
%% recorder(:,2) = node_size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf sq error
%% recorder(:,6) = mpf time
sizes = unique(K_data(:,1));

meanresults = nan(size(sizes, 1), 5);
stdresults = nan(size(sizes, 1), 5);

for i = 1:length(sizes)
    current_K_data = K_data(K_data(:,1) == sizes(i),:);
    K_mean_results(i, :) = mean(current_K_data);
    K_std_results(i, :) = std(current_K_data);
end


clf;
hold on;
grid on;
set(gca, 'xscale', 'log', 'yscale', 'log');
xlim([10^1.9, 10^4.1]);
ylim([10^-2.1, 10^3]);
axis square;

%% 4 x 4
errorbar(grid_mean_results(:, 1), grid_mean_results(:, 2), grid_std_results(:, 2), ...
         'o-', 'color', [0 0 0], 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         [0 0 0]);
errorbar(grid_mean_results(:, 1), grid_mean_results(:, 4), grid_std_results(:, 4), ...
         'o-', 'color', [0 0.5 0], 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         [0 0.5 0]);
%% Complete graph
errorbar(K_mean_results(:, 1), K_mean_results(:, 3), K_std_results(:, 3), ...
         'ro-', 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         'r');
errorbar(K_mean_results(:, 1), K_mean_results(:, 5), K_std_results(:, 5), ...
         'bo-', 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         'b');
xlabel('Sample size', 'interpreter', 'latex', 'fontsize', 25)
ylabel('Mean Squared Error', 'interpreter', 'latex', 'fontsize', 25)
h = legend('MLE Grid', 'MPF Grid', 'MLE K$_6$', 'MPF K$_6$', 'location', 'northeast');
set(h, 'interpreter', 'latex', 'fontsize', 25);

%% Save figure
saveas(gcf,'plots/effiency_4K','epsc')