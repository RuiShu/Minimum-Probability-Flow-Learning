function neighbor = is_neighbor(j, k, adj)
%% Check if nodes j and k are neighbors based on the connectivity matrix adj
%% 
%% Keywords
%% j   -- index of the first node
%% k   -- index of the second node
%% adj -- adjacency matrix

neighbor = adj(j, k) == 1