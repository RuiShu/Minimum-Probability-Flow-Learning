function main()

%% Load paths
addpath(genpath('./gen_utilities'));
addpath(genpath('./gen_utilities/minFunc_2012'));
addpath(genpath('./mle_utilities'));
addpath(genpath('./mpf_utilities'));

%% Settings for MSE computation
ntrials = 10;
edge_by_edge = 0;
mpf = 1;
mle = 1;
mpf_edge = 0;

%% Settings for data generation
grid_shape = [4 4];
nsamples = 10000;
adj = adj_mat_generation(grid_shape); % adjacency matrix

% burnin = 100*size(adj, 1);
% independent_steps = 10*size(adj, 1);
burnin = 1000;
independent_steps = 300;

%% Initialize data-recording variables
sq_err_Null = zeros(ntrials, 1);

if mpf
    sq_err_MPF = zeros(ntrials, 1);
    time_MPF = zeros(ntrials, 1);
end

if mpf_edge
    sq_err_MPF_edge = zeros(ntrials, 1);
    time_MPF_edge = zeros(ntrials, 1);
end

if mle
    sq_err_MLE = zeros(ntrials, 1);
    time_MLE = zeros(ntrials, 1);
end

J = ising_edge_generation(adj);

for trial = 1:ntrials % can use parfor loop for parallelization
    Xall_bin = ising_data_generation(J, nsamples, 'original', ...
                                     burnin, independent_steps); % for MPF

    fprintf(['Finished generating data Xall. size(Xall) = [# nodes, ' ...
             '# samples]:\n']);
    disp(size(Xall_bin));
    
    %% Minimum Probability Flow
    if mpf
        fprintf('MPF Estimation:\n');
        J_estimate = zeros(size(adj));
        
        t = tic();
        if edge_by_edge
            v = 1:size(J, 1);
            for j = 1:size(J, 1)
                for k = (j+1):size(J, 1) % v(adj(j, :) == 1)
                    J_estimate(j, k) = mpf_edge_estimation(Xall_bin, ...
                                                           adj, j, k, grid_shape);
                    J_estimate(k, j) = J_estimate(j, k);
                end
            end
        else
            J_estimate = mpf_estimate(adj, Xall_bin);
        end
        t = toc(t);
        time_MPF(trial) = t;
        
        fprintf('Square error: ');
        sq_err_MPF(trial) = sum((J(:) - J_estimate(:)).^2);
        sq_err_Null(trial) = sum(J(:).^2);
        disp(sq_err_MPF(trial));
    end

    %% MPF Edge
    if mpf_edge
        fprintf('MPF Edge Estimation:\n');
        J_estimate = zeros(size(adj));
        
        t = tic();
        v = 1:size(J, 1);
        for j = 1:size(J, 1)
            kv = v(adj(j, :) == 1);
            for k = kv(kv > j)
                J_estimate(j, k) = mpf_edge_estimation(Xall_bin, ...
                                                       adj, j, k, grid_shape);
                J_estimate(k, j) = J_estimate(j, k);
            end
        end
        t = toc(t);
        time_MPF_edge(trial) = t;
        
        fprintf('Square error: ');
        sq_err_MPF_edge(trial) = sum((J(:) - J_estimate(:)).^2);
        disp(sq_err_MPF_edge(trial));
    end

    if mle
        %% Maximum Likelihood Estimation
        fprintf('MLE Estimation:\n');
        t = tic();
        if edge_by_edge
            for j = 1:size(J, 1)
                for k = (j+1):size(J, 2)
                    J_estimate(j, k) = mle_edge_estimation(Xall_bin, adj, j, k);
                    J_estimate(k, j) = J_estimate(j, k);
                end
            end
        else
            J_estimate = mle_estimate(adj, Xall_bin);
        end
        t = toc(t);
        time_MLE(trial) = t;

        fprintf('Square error: ');
        sq_err_MLE(trial) = sum((J(:) - J_estimate(:)).^2);
        disp(sq_err_MLE(trial));
    end
end

%% Results
if mpf
    fprintf('MPF Average run time:      ');
    disp(mean(time_MPF))

    fprintf('MPF MSE:                   ');
    disp(mean(sq_err_MPF));
end

if mpf_edge
    fprintf('MPF Edge Average run time: ');
    disp(mean(time_MPF_edge))

    fprintf('MPF Edge MSE:              ');
    disp(mean(sq_err_MPF_edge));
end

if mle
    fprintf('MLE Average run time:      ');
    disp(mean(time_MLE))

    fprintf('MLE MSE:                   ');
    disp(mean(sq_err_MLE))
end 

fprintf('Null hypothesis MSE:       ');
disp(mean(sq_err_Null));

