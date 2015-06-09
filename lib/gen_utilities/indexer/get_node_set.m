function node_set = get_node_set(grid_shape, n_size, adj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the shape of the ising grid and the size of the "inner                          
% auxiliary MRF," determine the 2D indices of all corners of all inner aux MRF
% 
% Keywords:                                                                                
% grid_shape -- A (1,2) matrix describing the size of the lattice grid                     
% n_size     -- A numeric describing a single dimension of the inner aux MRF square grid.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assuming that the grid_shape is a square, and that grid_shape(1)
%% is a multiple of n_size, we can get the following corners:
idx = (0:(grid_shape(1)/n_size - 1))*n_size + 1;
idx = reshape(idx, length(idx), 1);

idx_2 = repmat(idx, [length(idx) 1]);
idx_1 = zeros(size(idx_2));

for i = 1:length(idx)
    idx_1(((i - 1)*length(idx) + 1):(i*length(idx))) = idx(i)*ones(length(idx), 1);
end

corners_idx_2d = [idx_1 idx_2];

%% Create a set of inner nodes and nodes
node_set = cell(length(idx), 2);

for i = 1:size(corners_idx_2d, 1)
    [inner_nodes, nodes] = get_nodes_by_corner(corners_idx_2d(i,:), ...
                                               grid_shape, n_size, adj);
    node_set{i, 1} = inner_nodes;
    node_set{i, 2} = nodes;
end
