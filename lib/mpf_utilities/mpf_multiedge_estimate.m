function Jnew_aux = mpf_multiedge_estimate(adj, Xall, inner_nodes, ...
                                           nodes) 
%% Find the aux indices of all inner nodes
v = 1:length(nodes);
aux_inner_idx = v(ismember(nodes, inner_nodes));

%% Extract auxiliary MRF. Use adjusted auxiliary adjacency method
d = length(nodes);
aux_adj = adj(nodes, nodes);
adjusted_aux_adj = ones(d) - eye(d);
adjusted_aux_adj(aux_inner_idx, :) = aux_adj(aux_inner_idx, :);
adjusted_aux_adj(:, aux_inner_idx) = aux_adj(:, aux_inner_idx);

%% Extract data
Xall_aux = Xall(nodes, :);

%% Construct adjusted aux Jnew
Jnew = randn(d, d) / sqrt(d) / 100;
Jnew = Jnew + Jnew';
Jnew_aux = Jnew .* adjusted_aux_adj;

%% minFunc Settings
% maxlinesearch is excessive just to be safe!!!!!! 
% learning works fine if this is just a few hundred 
maxlinesearch = 10000; 
% the number of gibbs sampling steps to take between samples

minf_options = [];
minf_options.display = 'none';
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;

%%%%%%%%%%% choose one of these two Ising model objective functions %%%%%%%%%%%%%
% K_dK_ising is slightly faster, and includes connectivity only to states
% which differ by a single bit flip.
%Jnew = minFunc( @K_dK_ising, Jnew(:), minf_options, Xall );
% K_dK_ising_allbitflipextension corresponds to a slightly modified choice
% of connectivity for MPF. This modified connectivity includes an
% additional connection between each state and the state which is reached
% by flipping all bits.  This connectivity pattern performs better in cases
% (like neural spike trains) where activity is extremely sparse.
Jnew_aux = minFunc( @K_dK_ising_allbitflipextension, Jnew_aux(:), minf_options, ...
                Xall_aux, adjusted_aux_adj);

Jnew_aux = reshape(Jnew_aux, size(adjusted_aux_adj)) .* adj(nodes, nodes);