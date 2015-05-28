function Jnew = mle_estimate(adj, Xall_bin)

d = size(adj, 1);
Jnew = randn( d, d ) / sqrt(d) / 100;
Jnew = Jnew + Jnew';
Jnew = Jnew/2 .* adj;

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
Jnew = minFunc( @L_dL_ising, Jnew(:), minf_options, Xall_bin, adj);

Jnew = reshape(Jnew, size(adj));