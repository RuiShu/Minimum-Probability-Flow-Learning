function ising_edge_estimation(j, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the edge weight between j and k in an ising lattice based on %
% multiple methods:                                                     %
% Full Ising model estimation method                                    %
% Clique extraction and estimation method                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize
nsamples = 1000;                 % number of training samples
adj = adj_mat_generation(4, 4); % adjacency matrix
[d, e] = size(adj);             % number of units

% Check if j and k are even connected in adj matrix
if adj(j, k) == 0
    fprintf('Not connected\n');
    return
end

% Find all nodes in the clique
aux_idx = (adj(j,:) + adj(k,:)) > 0;
% Index of j and k in the clique
aux_j = sum(aux_idx(1:j));
aux_k = sum(aux_idx(1:k));

%% minFunc Settings
addpath(genpath('./minFunc_2012'));
% maxlinesearch is excessive just to be safe!!!!!! 
% learning works fine if this is just a few hundred 
maxlinesearch = 10000; 
% the number of gibbs sampling steps to take between samples
independent_steps = 10*d; 

minf_options = [];
%options.display = 'none';
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;

run_checkgrad = 0;

if run_checkgrad
   d = 5;
   nsamples = 2;
   minf_options.DerivativeCheck = 'on';
else
   % make the weight matrix repeatable
%    rand('twister',355672);
%    randn('state',355672);
end

%% choose a random coupling matrix to generate the test data
J = randn( d, d ) / sqrt(d) * 3.;
J = J + J';
J = J/2;
J = J - diag(diag(J)); % set the diagonal so all the units are 0 bias
J = J - diag(sum(J));
J = adj .* J;

%% and generate the test data ...
fprintf( 'Generating %d training samples\n', nsamples );
burnin = 100*d;
t_samp = tic();
Xall = sample_ising( J, nsamples, burnin, independent_steps );
fprintf('Xall size \n');
disp(size(Xall))

%% generate auxiliary data set for 3 by 4 lattice grid...
Xall_aux = Xall(aux_idx,:);

t_samp = toc(t_samp);
fprintf( 'training sample generation in %f seconds \n', t_samp );

%% randomly initialize the parameter matrix we're going to try to learn
%% note that the bias units lie on the diagonal of J
Jnew = randn( d, d ) / sqrt(d) / 100;
Jnew = Jnew + Jnew';
Jnew = Jnew/2 .* adj;
Jnew_restricted = Jnew/2 .* adj;

% Construct adjacency matrix for auxiliary matrix
aux_adj = adj(aux_idx, aux_idx)
aux_adj_subrestricted = ones(size(aux_adj));
aux_adj_subrestricted([aux_j aux_k], :) = aux_adj([aux_j aux_k], :);
aux_adj_subrestricted(:, [aux_j aux_k]) = aux_adj(:, [aux_j aux_k]);
superaux_adj = [0 1 0; 1 0 1; 0 1 0];

% Construct J for auxiliary matrix
Jnew_aux = Jnew(aux_idx, aux_idx) .* aux_adj; % fully connected aux
Jnew_aux_restricted = Jnew_aux;         % restricted aux
Jnew_aux_subrestricted = Jnew_aux;
Jnew_superaux = Jnew([4 7 8], [4 7 8]) .* superaux_adj;

%% perform parameter estimation
fprintf( '\nRunning minFunc for up to %d learning steps...\n', maxlinesearch );
t_min = tic();

%%%%%%%%%%% choose one of these two Ising model objective functions %%%%%%%%%%%%%
% K_dK_ising is slightly faster, and includes connectivity only to states
% which differ by a single bit flip.
%Jnew = minFunc( @K_dK_ising, Jnew(:), minf_options, Xall );
% K_dK_ising_allbitflipextension corresponds to a slightly modified choice
% of connectivity for MPF. This modified connectivity includes an
% additional connection between each state and the state which is reached
% by flipping all bits.  This connectivity pattern performs better in cases
% (like neural spike trains) where activity is extremely sparse.

% Jnew = minFunc( @K_dK_ising_allbitflipextension, Jnew(:), minf_options, ...
%                 Xall, adj); 
% Jnew_aux = minFunc( @K_dK_ising_allbitflipextension, Jnew_aux(:), ...
%                     minf_options, Xall_aux, adj(aux_idx, aux_idx));
Jnew = minFunc( @K_dK_ising_allbitflipextension, Jnew(:), minf_options, ...
                Xall, ones(size(adj))); 
Jnew_restricted = minFunc( @K_dK_ising_allbitflipextension, ...
                           Jnew_restricted(:), minf_options, Xall, adj); 
Jnew_aux = minFunc( @K_dK_ising_allbitflipextension, Jnew_aux(:), ...
                    minf_options, Xall_aux, ones(size(aux_adj)));
Jnew_aux_restricted = minFunc( @K_dK_ising_allbitflipextension, Jnew_aux_restricted(:), ...
                    minf_options, Xall_aux, aux_adj);
Jnew_aux_subrestricted = minFunc( @K_dK_ising_allbitflipextension, ...
                                  Jnew_aux_subrestricted(:), minf_options, ...
                                  Xall_aux, aux_adj_subrestricted);
Jnew_superaux = minFunc( @K_dK_ising_allbitflipextension, ...
                                  Jnew_superaux(:), minf_options, ...
                                  Xall([4 7 8], :), superaux_adj);

Jnew = reshape(Jnew, size(J));
Jnew_restricted = reshape(Jnew_restricted, size(J));
Jnew_aux = reshape(Jnew_aux, sum(aux_idx), sum(aux_idx));
Jnew_aux_restricted = reshape(Jnew_aux_restricted, sum(aux_idx), sum(aux_idx));
Jnew_aux_subrestricted = reshape(Jnew_aux_subrestricted, sum(aux_idx), sum(aux_idx));
Jnew_superaux = reshape(Jnew_superaux, 3, 3);
t_min = toc(t_min);
fprintf( 'parameter estimation in %f seconds \n', t_min );

fprintf('True parameter                              ');
disp(J(j, k))
disp(J(5, 6))
disp(J(7, 8))
disp(J(2, 6))
disp(J(3, 7))
disp(J(6, 10))
disp(J(7, 11))
disp(J(2, 3))
disp(J(10, 11))
% fprintf('Estimated parameter by unrestricted full MRF');
% disp(Jnew(j, k))
% fprintf('Estimated parameter by restricted full MRF  ');
% disp(Jnew_restricted(j, k))
% fprintf('Estimated parameter by unrestricted aux MRF ');
% disp(Jnew_aux(aux_j, aux_k))
% fprintf('Estimated parameter by restricted aux MRF   ');
% disp(Jnew_aux_restricted(aux_j, aux_k))
fprintf('Estimated parameter by subrestricted aux MRF  ');
disp(Jnew_aux_subrestricted(aux_j, aux_k))
disp(Jnew_aux_subrestricted(3, 4))
disp(Jnew_aux_subrestricted(5, 6))
disp(Jnew_aux_subrestricted(1, 4))
disp(Jnew_aux_subrestricted(2, 5))
disp(Jnew_aux_subrestricted(4, 7))
disp(Jnew_aux_subrestricted(5, 8))
disp(Jnew_aux_subrestricted(1, 2))
disp(Jnew_aux_subrestricted(7, 8))
% fprintf('Estimated parameter by superaux MRF \n');
% disp(Jnew_superaux)



% fprintf('True parameter                              ');
% disp(J)
% fprintf('Estimated parameter by unrestricted full MRF');
% disp(Jnew)
% fprintf('Estimated parameter by restricted full MRF  ');
% disp(Jnew_restricted)
% fprintf('Estimated parameter by unrestricted aux MRF ');
% disp(Jnew_aux)
% fprintf('Estimated parameter by restricted aux MRF   ');
% disp(Jnew_aux_restricted)
% fprintf('Estimated parameter by subrestricted aux MRF');
% disp(Jnew_aux_subrestricted)
% fprintf('Estimated parameter by superaux MRF \n');
% disp(Jnew_superaux)

% fprintf('Control study')
% disp(J)
% disp(Jnew)
% disp(Jnew_aux)
% disp(Jnew_aux_restricted)