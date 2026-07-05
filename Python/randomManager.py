# randomManager.py

import numpy as np
from scipy.stats import beta, norm, lognorm, binom

class RandomNumberManager:
    def __init__(self, batch_size=1000000):
        self.batch_size = batch_size
        # One buffer & index per distribution type
        self._buffers = {}
        self._indices = {}

        for dist in [
            'beta', 'normal', 'lognormal', 'uniform', 'binomial', 'log10normal'
        ]:
            self._refill(dist)

    def _refill(self, dist):
        seg_edges = np.linspace(0, 1, self.batch_size + 1)
        buffer = np.random.uniform(low=seg_edges[:-1], high=seg_edges[1:])
        np.random.shuffle(buffer)
        self._buffers[dist] = buffer
        self._indices[dist] = 0

    def _get(self, dist):
        if self._indices[dist] >= self.batch_size:
            self._refill(dist)
        val = self._buffers[dist][self._indices[dist]]
        self._indices[dist] += 1
        return val

# Beta distribution with optional loc/scale for scaling (default [0,1])
def random_beta_lhs(rng, a, b, loc=0, scale=1):
    u = rng._get('beta')
    return beta.ppf(u, a, b) * scale + loc

def random_normal_lhs(rng, mu, sigma):
    u = rng._get('normal')
    return norm.ppf(u, loc=mu, scale=sigma)

def random_lognormal_lhs(rng, mean, sigma):
    u = rng._get('lognormal')
    return lognorm.ppf(u, s=sigma, scale=np.exp(mean))

def random_uniform_lhs(rng, low, high):
    u = rng._get('uniform')
    return low + (high - low) * u

def random_binomial_lhs(rng, n, p):
    u = rng._get('binomial')
    # Approximate inverse CDF for discrete binomial
    cdf = 0.0
    for k in range(n + 1):
        cdf += binom.pmf(k, n, p)
        if u <= cdf:
            return k
    return n

def random_log10normal_lhs(rng, mu, sigma):
    u = rng._get('log10normal')
    val = norm.ppf(u, loc=mu, scale=sigma)
    return 10 ** val
