function edge_weight = mle_edge_estimation(Xall, adj, j, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the edge weight between j and k in an ising lattice based on %
% multiple methods:                                                     %
% Full Ising model estimation method                                    %
% Clique extraction and estimation method                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize
d = size(adj, 1);                       % number of units

% Check if j and k are even connected in adj matrix
if adj(j, k) == 0
%     fprintf('Not connected\n');
    edge_weight = 0;
    return
end

% Find all nodes in the clique
aux_idx = (adj(j,:) + adj(k,:)) > 0;
% Index of j and k in the clique
aux_j = sum(aux_idx(1:j));
aux_k = sum(aux_idx(1:k));

%% minFunc Settings
% maxlinesearch is excessive just to be safe!!!!!! 
% learning works fine if this is just a few hundred 
maxlinesearch = 10000; 
% the number of gibbs sampling steps to take between samples

minf_options = [];
minf_options.display = 'none';
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;

% Note sure what this does
% run_checkgrad = 0;
% 
% if run_checkgrad
%     d = 5;
%     nsamples = 2;
%     minf_options.DerivativeCheck = 'on';
% else
% make the weight matrix repeatable
%    rand('twister',355672);
%    randn('state',355672);
% end

%% generate auxiliary data set for lattice grid...
Xall_aux = Xall(aux_idx,:);

%% randomly initialize the parameter matrix we're going to try to learn
%% note that the bias units lie on the diagonal of J
Jnew = randn( d, d ) / sqrt(d) / 100;
Jnew = Jnew + Jnew';
Jnew = Jnew/2 .* adj;                   % edge weights init
Jnew_restricted = Jnew/2 .* adj;        % impose lattice constraint

% Construct adjacency matrix for auxiliary matrix
aux_adj = adj(aux_idx, aux_idx);         % full aux matrix
aux_adj_subrestricted = ones(size(aux_adj)); 
aux_adj_subrestricted([aux_j aux_k], :) = aux_adj([aux_j aux_k], :);
aux_adj_subrestricted(:, [aux_j aux_k]) = aux_adj(:, [aux_j aux_k]); 

% Construct J for auxiliary matrix
Jnew_aux = Jnew(aux_idx, aux_idx);
Jnew_aux_restricted = Jnew_aux .* aux_adj;
Jnew_aux_subrestricted = Jnew_aux .* aux_adj_subrestricted;
% Jnew_aux = Jnew(aux_idx, aux_idx) .* aux_adj; % fully connected aux
% Jnew_aux_restricted = Jnew_aux;         % restricted aux
% Jnew_aux_subrestricted = Jnew_aux .* ;

%% perform parameter estimation
% fprintf( 'Running minFunc for up to %d learning steps...\n', maxlinesearch );
t_min = tic();

%%%%%%%%%%% choose one of these two Ising model objective functions %%%%%%%%%%%%%
% K_dK_ising is slightly faster, and includes connectivity only to states
% which differ by a single bit flip.
%Jnew = minFunc( @K_dK_ising, Jnew(:), minf_options, Xall );
% L_dL_ising corresponds to a slightly modified choice
% of connectivity for MPF. This modified connectivity includes an
% additional connection between each state and the state which is reached
% by flipping all bits.  This connectivity pattern performs better in cases
% (like neural spike trains) where activity is extremely sparse.

% method = 'Aux MRF Adjusted Restricted Connection';
method = 'AMAC';    
if strcmp(method, 'FMFC') % full MRF full connection
% MPF estimation over full Markov Random Field
    Jnew = minFunc( @L_dL_ising, Jnew(:), minf_options, ...
                    Xall, ones(size(adj))); 
    Jnew = reshape(Jnew, size(adj));
    edge_weight = Jnew(j, k);
elseif strcmp(method, 'FMRC') % full MRF restricted connection
% MPF estimation over full MRV with known lattice constraint
    Jnew_restricted = minFunc( @L_dL_ising, ...
                               Jnew_restricted(:), minf_options, Xall, ...
                               adj); 
    Jnew_restricted = reshape(Jnew_restricted, size(adj));
    edge_weight = Jnew_restricted(j, k);
elseif strcmp(method, 'AMFC') % aux MRF full connection
% MPF estimation over auxiliary MRV
    Jnew_aux = minFunc( @L_dL_ising, Jnew_aux(:), ...
                        minf_options, Xall_aux, ones(size(aux_adj)));
    Jnew_aux = reshape(Jnew_aux, sum(aux_idx), sum(aux_idx));
    edge_weight = Jnew_aux(aux_j, aux_k);
elseif strcmp(method, 'AMRC') % aux MRF restricted connection
% MPF estimation over auxiliary MRV with known lattice constraint
    Jnew_aux_restricted = minFunc( @L_dL_ising, ...
                                   Jnew_aux_restricted(:), ...
                                   minf_options, Xall_aux, aux_adj);
    Jnew_aux_restricted = reshape(Jnew_aux_restricted, sum(aux_idx), ...
                                  sum(aux_idx));
    edge_weight = Jnew_aux_restricted(aux_j, aux_k);
elseif strcmp(method, 'AMAC') % aux MRF adjusted connection
% MPF estimation over auxiliary MRV with adjusted lattice constraint (accounts for indirect connections in full MRF not present in auxiliary MRF)
    Jnew_aux_subrestricted = minFunc( @L_dL_ising, ...
                                      Jnew_aux_subrestricted(:), minf_options, ...
                                      Xall_aux, aux_adj_subrestricted);
    Jnew_aux_subrestricted = reshape(Jnew_aux_subrestricted, ...
                                     sum(aux_idx), sum(aux_idx));
    edge_weight = Jnew_aux_subrestricted(aux_j, aux_k);
end

% t_min = toc(t_min);
% fprintf( 'parameter estimation in %f seconds \n', t_min );

% Generate edgeweight from auxiliary subrestricted MRF

if false
    fprintf('True parameter                              \n'); 
    vertices = 1:(colsize*rowsize); 
    for i = vertices(adj(j,:) == 1)
        if i == k
            fprintf('>>>'); 
        else
            fprintf('   '); 
        end
        disp(J(j,i))
    end
    disp(J(j,k))

    fprintf('Estimated parameter by unrestricted full MRF\n'); 
    for i = vertices(adj(j,:) == 1)
        if i == k
            fprintf('>>>'); 
        else
            fprintf('   '); 
        end
        disp(Jnew(j,i))
    end
    disp(Jnew(j, k))

    fprintf('Estimated parameter by restricted full MRF  \n'); 
    for i = vertices(adj(j,:) == 1)
        if i == k
            fprintf('>>>'); 
        else
            fprintf('   '); 
        end
        disp(Jnew_restricted(j,i))
    end
    disp(Jnew_restricted(j, k))

    fprintf('Estimated parameter by unrestricted aux MRF \n'); 
    for i = vertices(aux_adj(aux_j,:) == 1)
        if i == aux_k
            fprintf('>>>'); 
        else
            fprintf('   '); 
        end
        disp(Jnew_aux(aux_j,i))
    end
    disp(Jnew_aux(aux_j, aux_k))

    fprintf('Estimated parameter by restricted aux MRF   \n'); 
    for i = vertices(aux_adj(aux_j,:) == 1)
        if i == aux_k
            fprintf('>>>'); 
        else
            fprintf('   '); 
        end
        disp(Jnew_aux_restricted(aux_j,i))
    end
    disp(Jnew_aux_restricted(aux_j, aux_k))

    fprintf('Estimated parameter by subrestricted aux MRF\n'); 
    for i = vertices(aux_adj(aux_j,:) == 1)
        if i == aux_k
            fprintf('>>>'); 
        else
            fprintf('   '); 
        end
        disp(Jnew_aux_subrestricted(aux_j,i))
    end
    disp(Jnew_aux_subrestricted(aux_j, aux_k))
end