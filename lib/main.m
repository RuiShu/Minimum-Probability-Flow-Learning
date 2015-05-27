function main()

%% Load paths
addpath(genpath('./gen_utilities'));
addpath(genpath('./gen_utilities/minFunc_2012'));
addpath(genpath('./mle_utilities'));
addpath(genpath('./mpf_utilities'));

%% Settings for MSE computation
ntrials = 1;
sq_err_Null = zeros(ntrials, 1);
sq_err_MPF = zeros(ntrials, 1);
sq_err_MLE = zeros(ntrials, 1);
time_MPF = zeros(ntrials, 1);
time_MLE = zeros(ntrials, 1);

%% Settings for data generation
grid_shape = [3 1];
nsamples = 10000;
adj = adj_mat_generation(grid_shape); % adjacency matrix

%% Adjustment
adj = ones(3) - eye(3);
J = ising_edge_generation(adj);
%% J = [0 -1.4 ; -1.4 0];
disp(J);

for trial = 1:ntrials % can use parfor loop for parallelization
%     Xall_bin = (Xall_dic + 1)/2;
%     Xall_bin = ising_data_generation(J/2, nsamples);
%     Xall_dic = (Xall_bin - 0.5)*2;
%     disp(mean(Xall_dic, 2))
%     disp(mean(Xall_dic(1, :) .* Xall_dic(2, :), 2))

%     Xall_dic1 = ising_data_generation_james(J, nsamples);
%     Xall_bin1 = (Xall_dic + 1)/2;
%     disp(mean(Xall_dic1, 2));
%     disp(mean(Xall_dic1(1, :) .* Xall_dic1(2, :), 2))

    Xall_dic = ising_data_generation(J, nsamples, 'new'); % for James
    Xall_bin = ising_data_generation(J, nsamples, 'original'); % for MPF
    fprintf('Info on Xall, direct dichotomous\n');
    fprintf('Info on Xall, binary-based dichotomous\n');
    Xall_dic2 = (Xall_bin - 0.5)*2;
    disp(Xall_dic2(:,1:10))
    disp(mean(Xall_dic2, 2))
    disp(mean(Xall_dic2(1, :) .* Xall_dic2(2, :), 2))

    fprintf(['Finished generating data Xall. size(Xall) = [# nodes, ' ...
             '# samples]:\n']);
    disp(size(Xall_bin));
    
    %% Minimum Probability Flow
    fprintf('MPF Estimation:\n');
    J_estimate = zeros(size(adj));

    t = tic();

    if 1
        J_estimate = mpf_estimate(adj, Xall_bin);
    else
        for j = 1:size(J, 1)
            for k = (j+1):size(J, 2)
                J_estimate(j, k) = ising_edge_estimation(Xall_bin, adj, j, k);
                J_estimate(k, j) = J_estimate(j, k);
            end
            fprintf('%.3f complete. Trial %3d out of %3d\n', ...
                    j/size(J, 1), trial, ntrials);
        end
    end

    t = toc(t);
    time_MPF(trial) = t;
    sq_err_MPF(trial) = sum((J(:) - J_estimate(:)).^2);
    sq_err_Null(trial) = sum(J(:).^2);

    disp(J_estimate);

    fprintf('Square error: ');
    disp(sum((J(:) - J_estimate(:)).^2));

    %% Maximum Likelihood Estimation
    fprintf('MLE Estimation:\n');
    t = tic();
    [theta_est, e] = variational_estimate(J, Xall_bin', ...
                                          1, 1:size(J, 1));
    t = toc(t);
    time_MLE(trial) = t;
    sq_err_MLE(trial) = e;

    
    fprintf('Square error: ');
    disp(e);
end

%% Results
fprintf('MPF Average run time:');
disp(mean(time_MPF))

fprintf('MPF MSE:             ');
disp(mean(sq_err_MPF));

fprintf('MLE Average run time:');
disp(mean(time_MLE))

fprintf('MLE MSE:             ');
disp(mean(sq_err_MLE))

fprintf('Null hypothesis MSE: ');
disp(mean(sq_err_Null));

