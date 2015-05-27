function Xall = ising_data_generation(J, d, nsamples)

%% and generate the test data ...
independent_steps = 10*d; 
fprintf( 'Generating %d training samples\r', nsamples);
burnin = 100*d;
t_samp = tic();
Xall = sample_ising( J, nsamples, burnin, independent_steps);
