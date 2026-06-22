"""mcmc_toy_model.py — 2-param MH-MCMC demo (beta0, z_star)"""
import numpy as np
import pandas as pd
import os

def logistic_pdf(z, mu, scale):
    exp_term = np.exp(-(z - mu) / scale)
    return exp_term / (scale * (1 + exp_term) ** 2)

def synthetic_data(beta0_true, z_star_true, sigma_true, noise_std, num_points=20):
    z_obs = np.linspace(0.0, 1.0, num_points)
    w_true = -1.0 - beta0_true * logistic_pdf(z_obs, z_star_true, sigma_true)
    return z_obs, w_true + np.random.normal(0.0, noise_std, size=num_points)

def log_likelihood(z, w_obs, beta0, z_star, sigma, noise_std):
    w_model = -1.0 - beta0 * logistic_pdf(z, z_star, sigma)
    return -0.5 * np.sum(((w_obs - w_model) / noise_std) ** 2)

def metropolis_hastings(z_obs, w_obs, sigma_true, noise_std, n_steps=20000, burn_in=5000,
                        init_params=(0.1, 0.3), proposal_scales=(0.01, 0.02)):
    beta0, z_star = init_params
    chain = np.zeros((n_steps, 2))
    current = log_likelihood(z_obs, w_obs, beta0, z_star, sigma_true, noise_std)
    for i in range(n_steps):
        b, zs = np.random.normal(beta0, proposal_scales[0]), np.random.normal(z_star, proposal_scales[1])
        if b <= 0 or zs < 0 or zs > 1:
            chain[i] = [beta0, z_star]; continue
        new = log_likelihood(z_obs, w_obs, b, zs, sigma_true, noise_std)
        if np.log(np.random.rand()) < (new - current):
            beta0, z_star, current = b, zs, new
        chain[i] = [beta0, z_star]
    return pd.DataFrame(chain[burn_in:], columns=['beta0', 'z_star'])

def main():
    np.random.seed(0)
    z_obs, w_obs = synthetic_data(0.1, 0.4, 0.2, 0.005, num_points=15)
    chain_df = metropolis_hastings(z_obs, w_obs, 0.2, 0.005, n_steps=12000, burn_in=2000,
                                   init_params=(0.05, 0.5), proposal_scales=(0.005, 0.01))
    summary = chain_df.describe().loc[['mean', 'std']]
    summary.to_csv('mcmc_summary.csv')
    print(summary)

if __name__ == '__main__':
    main()
