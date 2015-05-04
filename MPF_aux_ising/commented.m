fprintf( '\nGenerating samples using learned parameters for comparison...\n' );
Xnew = sample_ising( Jnew, nsamples, burnin, independent_steps );
fprintf( 'sample generation with learned parameters in %f seconds \n', t_samp );

% generate correlation matrices for the original and recovered coupling matrices
mns = mean( Xall, 2 );
Xt = Xall - mns(:, ones(1,nsamples));
sds = sqrt(mean( Xt.^2, 2 ));
Xt = Xt./sds(:, ones(1,nsamples));
Corr = Xt*Xt'/nsamples;
mns = mean( Xnew, 2 );
Xt = Xnew - mns(:, ones(1,nsamples));
sds = sqrt(mean( Xt.^2, 2 ));
Xt = Xt./sds(:, ones(1,nsamples));
Corrnew = Xt*Xt'/nsamples;

Jdiff = J - Jnew;
Corrdiff = Corr - Corrnew;
jmn = min( [J(:); Jnew(:); Jdiff(:)] );
jmx = max( [J(:); Jnew(:); Jdiff(:)] );
cmn = min( [Corr(:); Corrnew(:); Corrdiff(:)] );
cmx = max( [Corr(:); Corrnew(:); Corrdiff(:)] );

% show the original, recovered and differences in coupling matrices
figure();
subplot(2,3,1);
imagesc( J, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{ }' );
subplot(2,3,2);
imagesc( Jnew, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{new}' );
subplot(2,3,3);
imagesc( Jdiff, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{ } - {J}_{new}' );

% show the original, recovered and differences in correlation matrices
subplot(2,3,4);
imagesc( Corr, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{ }' );    
subplot(2,3,5);
imagesc( Corrnew, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{new}' );    
subplot(2,3,6);
imagesc( Corrdiff, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{ } - {C}_{new}' );    

figure();
plot( Corr(:), Corrnew(:), '.' );
axis([cmn,cmx,cmn,cmx]);
axis image;
xlabel( '{C}_{ }' );
ylabel( '{C}_{new}' );
title( 'scatter plot of correlations for original and recovered parameters' );
