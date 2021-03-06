function adj = adj_mat_generation(grid_shape)

r = grid_shape(1);
c = grid_shape(2);

diagVec1 = repmat([ones(c-1,1); 0],r,1);  %# Make the first diagonal vector
                                          %#   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);             %# Remove the last value
diagVec2 = ones(c*(r-1),1);               %# Make the second diagonal vector
                                          %#   (for vertical connections)
adj = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
      diag(diagVec2,c);
adj = adj+adj.';                         %'# Add the matrix to a transposed
                                          %#   copy of itself to make it
                                          %#   symmetric

