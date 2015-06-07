function J = ising_edge_generation(adj)

[d, ~] = size(adj);             % number of units
J = randn( d, d ) / sqrt(d) * 1.;
J = J + J';
J = J/2;
J = J - diag(diag(J)); % set the diagonal so all the units are 0 bias
J = J - diag(sum(J));
J = adj .* J;


