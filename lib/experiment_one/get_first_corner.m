function idx_2d = get_first_corner(grid_shape, n_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the shape of the ising grid and the size of the "inner                         %
% auxiliary MRF," determine the 2D index of the top left most corner                      %
% of the inner auxiliary MRF.                                                             %
% Definitions of inner auxiliary MRF: we a given auxiliary matrix, the nodes              %
% of interest are flanked by the neighborhood, which is necessary for proper estimation   %
% of the edges between the nodes of interest. the inner aux MRF is just a square grid.    %
%                                                                                         %
% Keywords:                                                                               %
% grid_shape -- A (1,2) matrix describing the size of the lattice grid                    %
% n_size     -- A numeric describing a single dimension of the inner aux MRF square grid. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_1 = floor((grid_shape(1) - n_size)/2);
idx_2 = floor((grid_shape(2) - n_size)/2);

idx_2d = [idx_1 idx_2]
