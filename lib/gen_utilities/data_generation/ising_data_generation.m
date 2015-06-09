function Xall = ising_data_generation(J, nsamples, method, ...
                                      burnin, independent_steps)

%% and generate the test data ...
fprintf( 'Generating %d training samples\r', nsamples);

if strcmp(method, 'original')
    Xall = sample_ising_original( J, nsamples, burnin, independent_steps);    
elseif strcmp(method, 'new')
    Xall = sample_ising( J, nsamples, burnin, independent_steps);    
end
