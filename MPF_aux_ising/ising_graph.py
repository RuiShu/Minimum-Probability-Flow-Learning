import numpy as np
from numpy import atleast_2d as vec

class IsingGraph(object):

    def __init__(self, W):
        """Initialization of IsingGraph object
        
        Keyword arguments:
        W -- edge weight matrix for ising model
        """
        self.__W = W

    def gibbs_sampler(self, n_samples, burn_in, independent_steps):
        W = self.__W 
        n_sampling_steps = burn_in + (n_samples - 1)*independent_steps
        n_dims = W.shape[0]

        # choose dimensions to gibbs sample
        upd_i = np.floor(vec(np.random.uniform(size=n_sampling_steps)).T*n_dims)

        # Precalculate the random numbers for comparison
        uni_rand = np.random.uniform(n_sampling_steps)
        X_out = np.zeros(shape=(n_dims, 1))

        i_out = 1
        next_sample = burn_in

        bias = np.diag(W)
        W = W - np.diag(bias)

        for si in range(0, n_sampling_steps):
            E_act = 2*W[upd_i[si,0],:]



if __name__ == "__main__":

    W = np.array([[1, 1], [1, 1]])
    graph = IsingGraph(W)

    graph.gibbs_sampler(3, 3, 3)
        
