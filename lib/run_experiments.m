function run_experiments()

%% Load paths
addpath(genpath('./gen_utilities'));
addpath(genpath('./gen_utilities/minFunc_2012'));
addpath(genpath('./mle_utilities'));
addpath(genpath('./mpf_utilities'));
addpath(genpath('./experiment_one'));

%% Run experiments
experiment_one()
