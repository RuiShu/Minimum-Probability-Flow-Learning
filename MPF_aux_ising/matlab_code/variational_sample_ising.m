function X = variational_sample_ising(T, niters, nsamples)
X = zeros(nsamples, size(T, 1));
for i = 1:nsamples
    X(i, :) = individual_sample(T, niters);
end


function x = individual_sample(T, niters)
d = size(T,1);
x = floor(2 * rand(1,d));x(x == 0) = -1;
for i = 1:niters
    index = randperm(d, 1);
    xp = x;xp(index) = xp(index) * -1;
    prob_x = unnormalized_probability(x, T);
    prob_xp = unnormalized_probability(xp, T);
    
    prob_transition = min(1, prob_xp / prob_x);
    if rand < prob_transition
        x = xp;
    end
end
