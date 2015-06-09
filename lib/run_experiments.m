function run_experiments()

%% Load paths
addpath(genpath('./gen_utilities/data_generation'));
addpath(genpath('./gen_utilities/indexer'));
addpath(genpath('./gen_utilities/minFunc_2012'));
addpath(genpath('./mle_utilities'));
addpath(genpath('./mpf_utilities'));
addpath(genpath('./experiment_1'));
addpath(genpath('./experiment_2'));
addpath(genpath('./experiment_3'));
addpath(genpath('./miniexperiment_1'));
addpath(genpath('./miniexperiment_2b'));

run_experiment_1 = 0;
run_experiment_2 = 0;
run_experiment_3 = 1;
run_miniexperiment_1 = 0;
run_miniexperiment_2a = 0;
run_miniexperiment_2b = 0;

%% Note: remember to change how J is generated?

%% Run experiments

%% experiment_1: For 4x4 grid, for n = [10 100 1000 10000], run
%% mpf v. mle.
%% Graph type: ising grid, 4^2
%% Algorithm type: vanilla mle, mpf
%% Sample size: [10 100 1000 10000]
%% Comparing: algorithm speed, estimator efficiency
if run_experiment_1
    n = [10 100 1000 10000];
    add = [400 600 4000 6000]
    for nsamples = add
        experiment_1(nsamples)
    end
end

%% miniexperiment_1: For n = [10000], see how big of a sq grid mle
%% and mpf can analyze.
%% Graph type: ising grid, [1:10]^2
%% Algorithm type: vanilla mle, mpf
%% Sample size: 1000
%% Comparing: tractability
if run_miniexperiment_1
    dims = 1:5;
    for dim = dims
        miniexperiment_1(dim)
    end
end

%% miniexperiment_2a: For n = [10000], check if multi-patch node
%% selection yields good estimation, so long as neighborhood is
%% chosen correctly.
%% Graph type: ising grid, 10x10
%% Algorithm type: factorized mpf
%% Sample size: 10000
%% Comparing: theory validation

%% miniexperiment_2b: For n = [10000], compare various grid sizes
%% for mle and mpf. And find the optimal error/time trade off.
%% Graph type: ising grid, 30x30
%% Algorithm type: factorized mle, mpf
%% Factorization size: [1 2 3 5 6 10 15 30]
%% Sample size: 10000
%% Comparing: optimization

if run_miniexperiment_2b
    grid_shape = [30 30];
    n = [1 2 3 5 6 10 15];

    %% Check that grid_shape is square
    if grid_shape(1) ~= grid_shape(2)
        fprintf('Not a square\n');
        return
    end
    
    %% Check that all values in n are factors of grid_shape dim
    if any(bsxfun(@mod, grid_shape(1), n) > 0)
        fprintf('Not a value n_size partitioning\n');
        return
    end
    
    for n_size = n
        miniexperiment_2b(n_size)
    end
end

%% experiment_two: For 30x30 grid, for n = [10 100 1000 10000],
%% compare mpf v. mle. Same type of grid? Different type of grid?
%% Use different.

if run_experiment_2
    n = [10 100 1000 10000];

    for nsamples = n
        experiment_2(nsamples)
    end
end


%% experiment_three: complete graph, for n = [10 100 1000 10000],
%% run mpf v. mle for as large a graph as possible
if run_experiment_3
    n = [10 100];
    n = [1000 10000 13000 16000 19000 22000];
    n = [10000 13000 16000 19000 22000];
    nodes = 4:19;
    for nsamples = n
        for node_size = nodes
            experiment_3(nsamples, node_size)
        end
    end
end


