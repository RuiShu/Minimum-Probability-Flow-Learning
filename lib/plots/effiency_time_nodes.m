function effiency_time_nodes()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get and clean the data for both miniexperiment 2b n1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_file = ['miniexperiment_2b_factorization' int2str(1) '.mat'];
data = load(load_file);
data = data.recorder;

for i = [2 3 6 10 15 30]
    load_file = ['miniexperiment_2b_factorization' int2str(i) '.mat'];
    current_data = load(load_file);
    data = [data ; current_data.recorder];
end

%%%%%%%%%%%%%%%%%%%%%%%
% Construct the plots %
%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize data-recording variables
%% recorder is (ntrials, 6) matrix:
%% recorder(:,1) = sample size
%% recorder(:,2) = factorization size
%% recorder(:,3) = mle sq error
%% recorder(:,4) = mle time
%% recorder(:,5) = mpf sq error
%% recorder(:,6) = mpf time
node_sizes = unique(data(:, 2));

mean_results = nan(size(node_sizes, 1), size(data, 2));
std_results = nan(size(node_sizes, 1), size(data, 2));

for i = 1:length(node_sizes)
    current_data = data(data(:,2) == node_sizes(i),:);
    mean_results(i, :) = mean(current_data);
    std_results(i, :) = std(current_data);
end

%% Fitting a curve for MPF
X = [ones(size(data, 1), 1) data(:, 5) data(:, 5).^2 data(:, 5).^3 ...
     ];
Y = data(:, 6);
beta = X \ Y;
x_mpf = 11:0.01:14.1;
X_pred = [ones(length(x_mpf), 1) x_mpf' x_mpf'.^2 x_mpf'.^3];
y_mpf = X_pred * beta;

%% Fitting a curve for MLE
idx = ~isnan(data(:, 3));
X = [ones(sum(idx), 1) data(idx, 3)];
Y = data(idx, 4);
beta = X \ Y;
x_mle = 11:0.01:14.1;
X_pred = [ones(length(x_mle), 1) x_mle'];
y_mle = X_pred * beta;

clf;
hold on;
grid on;
xlim([10.7 14.5]);
ylim([3 10]);
axis square;

plot(x_mle, y_mle, 'r', 'linewidth', 2);
plot(x_mpf, y_mpf, 'b', 'linewidth', 2);
plot(data(:, 3), data(:, 4), 'rx', 'markersize', 10, 'linewidth', 2); 
plot(data(:, 5), data(:, 6), 'bx', 'markersize', 10, 'linewidth', 2);
xlabel('Mean Squared Error', 'interpreter', 'latex', 'fontsize', 25)
ylabel('Time (s)', 'interpreter', 'latex', 'fontsize', 25)
h = legend('MLE Grid', 'MPF Grid', 'location', 'northeast');
set(h, 'interpreter', 'latex', 'fontsize', 25);

%% Save figure
saveas(gcf,'plots/effiency_time_nodes','epsc')