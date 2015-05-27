function Xall = ising_data_generation_james(J, nsamples)

%% and generate the test data ...
fprintf( 'Generating %d training samples\r', nsamples);
niters = 100*size(J, 1);
Xall = variational_sample_ising(J, niters, nsamples)';
