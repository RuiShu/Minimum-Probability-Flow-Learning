function time_4()
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get and clean the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4 by 4 grid
load_file = ['experiment_1_n' int2str(10) '.mat'];
grid_data = load(load_file);
grid_data = grid_data.recorder;

for i = [100 1000 10000]
    load_file = ['experiment_1_n' int2str(i) '.mat'];
    current_grid_data = load(load_file);
    grid_data = [grid_data ; current_grid_data.recorder];
end

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

clf;
hold on;
grid on;
set(gca, 'xscale', 'log', 'yscale', 'log');
xlim([0, 10^4.1]);
ylim([10^-2.2, 10^.4]);
axis square

%% 4 x 4
errorbar(grid_mean_results(:, 1), grid_mean_results(:, 3), grid_std_results(:, 3), ...
         'ro--', 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         'r');
errorbar(grid_mean_results(:, 1), grid_mean_results(:, 5), grid_std_results(:, 5), ...
         'bo--', 'linewidth', 1, 'markersize', 5, 'markerfacecolor', ...
         'b');
xlabel('Sample size', 'interpreter', 'latex', 'fontsize', 25)
ylabel('Time (s)', 'interpreter', 'latex', 'fontsize', 25)
h = legend('MLE Grid', 'MPF Grid', 'location', 'northeast');
set(h, 'interpreter', 'latex', 'fontsize', 25);

%% Save figure
saveas(gcf,'plots/time_4','epsc')