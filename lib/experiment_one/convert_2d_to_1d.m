function nodes = convert_2d_to_1d(nodes_2d, grid_size)

nodes = nodes_2d(:, 1)*grid_size(1) + nodes_2d(:, 2);