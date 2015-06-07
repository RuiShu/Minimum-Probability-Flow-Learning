
function [inner_nodes, nodes] = get_all_nodes(grid_shape, n_size, adj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the shape of the ising grid and the size of the "inner                          
% auxiliary MRF," determine the 2D indices of all the nodes, both                          
% within the inner aux MRF as well as its neighborhood.  Definitions                       
% of inner auxiliary MRF: we a given auxiliary matrix, the nodes of                        
% interest are flanked by the neighborhood, which is necessary for                         
% proper estimation of the edges between the nodes of interest. the                        
% inner aux MRF is just a square grid.                                                     
% 
% Keywords:                                                                                
% grid_shape -- A (1,2) matrix describing the size of the lattice grid                     
% n_size     -- A numeric describing a single dimension of the inner aux MRF square grid.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine first corner
idx_1 = floor((grid_shape(1) - n_size)/2);
idx_2 = floor((grid_shape(2) - n_size)/2);
corner_idx_2d = [idx_1 idx_2];

%% Figure out all offsets
offset_2 = repmat([0:(n_size - 1)]', [n_size, 1]);
offset_1 = zeros(size(offset_2));

for i = 0:(n_size - 1)
    idx = i*n_size;
    offset_1((idx+1):(idx+n_size)) = i*ones(n_size, 1);
end

offset = [offset_1 offset_2];

%% Figure out all inner nodes
inner_nodes_2d = bsxfun(@plus, offset, corner_idx_2d);
inner_nodes = convert_2d_to_1d(inner_nodes_2d, grid_shape);

%% Figure out all nodes
nodes = (1:size(adj, 1));
nodes = nodes(sum(adj(inner_nodes,:), 1) > 0)';
