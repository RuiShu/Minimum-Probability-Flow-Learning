function Xall = ising_data_generation(J, nsamples, method)

%% and generate the test data ...
d = size(J, 1);
independent_steps = 10*d; 
fprintf( 'Generating %d training samples\r', nsamples);
burnin = 100*d;

if strcmp(method, 'original')
    Xall = sample_ising_original( J, nsamples, burnin, independent_steps);    
elseif strcmp(method, 'new')
    Xall = sample_ising( J, nsamples, burnin, independent_steps);    
end
