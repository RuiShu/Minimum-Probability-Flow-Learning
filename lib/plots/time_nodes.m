function time_nodes()

%% K_n time plots, samples = 10000
load_file = ['experiment_3_n10000_nodes' int2str(9) '.mat'];
data = load(load_file);
data = data.recorder;

for i = 10:19
    load_file = ['experiment_3_n10000_nodes' int2str(i) '.mat'];
    current_data = load(load_file);
    data = [data ; current_data.recorder];
end

%% Initialize data-recording variables
%% recorder is (ntrials, 6) matrix:
%% recorder(:,1) = nsamples
%% recorder(:,2) = node_size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf sq error
%% recorder(:,6) = mpf time
node_sizes = unique(data(:,2));

meanresults = nan(size(node_sizes, 1), size(data, 2));
stdresults = nan(size(node_sizes, 1), size(data, 2));

for i = 1:length(node_sizes)
    current_data = data(data(:,2) == node_sizes(i),:);
    mean_results(i, :) = mean(current_data);
    std_results(i, :) = std(current_data);
end

clf;
hold on;
grid on;
axis square;

mean_results = mean_results(1:end, :);
std_results = std_results(1:end, :);

%% 4 x 4
errorbar(mean_results(:, 2), mean_results(:, 4), std_results(:, 4), ...
         'ro--', 'linewidth', 1, 'markersize', 3, 'markerfacecolor', ...
         'r');
errorbar(mean_results(:, 2), mean_results(:, 6), std_results(:, 6), ...
         'bo--', 'linewidth', 1, 'markersize', 3, 'markerfacecolor', ...
         'b');
xlabel('Number of nodes (n)', 'interpreter', 'latex', 'fontsize', 25)
ylabel('Time (s)', 'interpreter', 'latex', 'fontsize', 25)
h = legend('MLE K$_n$', 'MPF K$_n$', 'location', 'northeast');
set(h, 'interpreter', 'latex', 'fontsize', 25);

%% Save figure
saveas(gcf,'plots/time_nodes','epsc')