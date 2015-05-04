function ising_edge_estimation(j, k)
% Computes the edge weight between j and k in an ising lattice based on
% multiple methods:
% Full Ising model estimation method
% Clique extraction and estimation method

addpath(genpath('./minFunc_2012'));

% initialize
adj = adj_mat_generation(3, 4);
[d, e] = size(adj); % number of units

nsamples = 100000; % number of training samples

maxlinesearch = 10000; % this number is excessive just to be safe!!!!!! learning works fine if this is just a few hundred
independent_steps = 10*d; % the number of gibbs sampling steps to take between samples

run_checkgrad = 0;

minf_options = [];
%options.display = 'none';
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;

if run_checkgrad
   d = 5;
   nsamples = 2;
   minf_options.DerivativeCheck = 'on';
else
   % make the weight matrix repeatable
   rand('twister',355672);
   randn('state',355672);
end

% choose a random coupling matrix to generate the test data
J = randn( d, d ) / sqrt(d) * 3.;
J = J + J';
J = J/2;
J = J - diag(diag(J)); % set the diagonal so all the units are 0 bias
J = J - diag(sum(J));
J = adj .* J;
fprintf( 'Generating %d training samples\n', nsamples );
burnin = 100*d;

% and generate the test data ...
t_samp = tic();
Xall = sample_ising( J, nsamples, burnin, independent_steps );
fprintf('Xall size \n');
disp(size(Xall))

% generate corrupt data set for 3 by 4 lattice grid...
Xall_corr = Xall;
Xall_corr(1,:) = round(rand(size(Xall(1, :))));
Xall_corr(4,:) = round(rand(size(Xall(1, :))));
Xall_corr(9,:) = round(rand(size(Xall(1, :))));
Xall_corr(12,:) = round(rand(size(Xall(1, :))));

Xall_aux = Xall([2,3,5,6,7,8,10,11],:);

t_samp = toc(t_samp);
fprintf( 'training sample generation in %f seconds \n', t_samp );

% randomly initialize the parameter matrix we're going to try to learn
% note that the bias units lie on the diagonal of J
Jnew = randn( d, d ) / sqrt(d) / 100;
Jnew = Jnew + Jnew';
Jnew = Jnew/2;
Jnew_corr = Jnew;
Jnew_aux = Jnew([2,3,5,6,7,8,10,11], [2,3,5,6,7,8,10,11]);

% perform parameter estimation
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

Jnew = minFunc( @K_dK_ising_allbitflipextension, Jnew(:), minf_options, Xall ); 
Jnew_corr = minFunc( @K_dK_ising_allbitflipextension, Jnew_corr(:), ...
                     minf_options, Xall_corr );

Jnew_aux = minFunc( @K_dK_ising_allbitflipextension, Jnew_aux(:), minf_options, Xall_aux);



Jnew = reshape(Jnew, size(J));
Jnew_corr = reshape(Jnew_corr, size(J));
Jnew_aux = reshape(Jnew_aux, size(J)-4);
t_min = toc(t_min);
fprintf( 'parameter estimation in %f seconds \n', t_min );

disp(J(6,7))
disp(Jnew(6,7))
disp(Jnew_corr(6,7))
disp(Jnew_aux(4,5))

% fprintf('Control study')
% disp(J)
% disp(Jnew)
% disp(Jnew_corr)
% disp(Jnew_aux)