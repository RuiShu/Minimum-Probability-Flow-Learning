function [inner_nodes, nodes] = get_nodes_by_index(index, n_size, adj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Keywords:                                                                                
% grid_shape -- A (1,2) matrix describing the size of the lattice grid                     
% n_size     -- A numeric describing a single dimension of the inner aux MRF square grid.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inner_nodes = [index:(index + n_size - 1)]';

%% Figure out all nodes
biased_adj = adj + eye(size(adj));
nodes = (1:size(adj, 1));
nodes = nodes(sum(biased_adj(inner_nodes,:), 1) > 0)';
