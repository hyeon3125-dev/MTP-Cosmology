"""advanced_mcmc_analysis.py — 3-param Metropolis-Hastings MCMC (beta0, z_star, sigma)"""
from __future__ import annotations
import numpy as np
import pandas as pd
import os

def logistic_pdf(z, mu, scale):
    if scale <= 0:
        raise ValueError("scale must be positive for logistic PDF")
    exp_term = np.exp(-(z - mu) / scale)
    return exp_term / (scale * (1 + exp_term) ** 2)

def generate_synthetic_data(beta0_true, z_star_true, sigma_true, z_min, z_max, num_points, noise_std):
    z_obs = np.linspace(z_min, z_max, num_points)
    w_true = -1.0 - beta0_true * logistic_pdf(z_obs, z_star_true, sigma_true)
    return z_obs, w_true + np.random.normal(0.0, noise_std, size=num_points)

def log_likelihood(z, w_obs, beta0, z_star, sigma, noise_std):
    if beta0 <= 0 or sigma <= 0: return -np.inf
    if z_star < 0 or z_star > np.max(z): return -np.inf
    w_model = -1.0 - beta0 * logistic_pdf(z, z_star, sigma)
    return -0.5 * np.sum(((w_obs - w_model) / noise_std) ** 2)

def run_mcmc(z_obs, w_obs, noise_std, n_steps, burn_in, init_params, proposal_scales):
    beta0, z_star, sigma = init_params
    chain = np.zeros((n_steps, 3))
    z_max = np.max(z_obs)
    current = log_likelihood(z_obs, w_obs, beta0, z_star, sigma, noise_std)
    for i in range(n_steps):
        b, zs, sg = (np.random.normal(beta0, proposal_scales[0]),
                     np.random.normal(z_star, proposal_scales[1]),
                     np.random.normal(sigma, proposal_scales[2]))
        if b <= 0 or sg <= 0 or zs < 0 or zs > z_max:
            chain[i] = [beta0, z_star, sigma]; continue
        new = log_likelihood(z_obs, w_obs, b, zs, sg, noise_std)
        if np.log(np.random.rand()) < (new - current):
            beta0, z_star, sigma, current = b, zs, sg, new
        chain[i] = [beta0, z_star, sigma]
    return pd.DataFrame(chain[burn_in:], columns=['beta0', 'z_star', 'sigma'])

def main():
    np.random.seed(42)
    z_obs, w_obs = generate_synthetic_data(0.1, 0.4, 0.3, 0.0, 1.5, 25, 0.01)
    chain_df = run_mcmc(z_obs, w_obs, 0.01, 20000, 4000, (0.05, 0.6, 0.2), (0.005, 0.02, 0.01))
    summary = chain_df.describe().loc[['mean', 'std']]
    chain_df.to_csv('advanced_mcmc_chain.csv', index=False)
    summary.to_csv('advanced_mcmc_summary.csv')
    print(summary)

if __name__ == '__main__':
    main()
