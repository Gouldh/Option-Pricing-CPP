#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>

namespace br = boost::random;

// -----------------------------------------------------------------------------------
// Author: Hunter Gould
// Date: 11/26/2023
// Description: This project is an experiment in translating past projects from Python
//              to C++. It involves the implementation of various option pricing models,
//              including the Black-Scholes Model, the Heston Model, and the Merton Jump
//              Diffusion Model. This translation is part of an effort to explore the
//              nuances and performance differences between Python and C++ in the context
//              of financial modeling and computational finance.
// -----------------------------------------------------------------------------------

// Normal CDF
double normCDF(double value) {
    return 0.5 * std::erfc(-value * std::sqrt(0.5));
}

// Black-Scholes Model
double blackScholes(double S, double X, double T, double r, double sigma) {
    double d1 = (std::log(S / X) + (r + 0.5 * std::pow(sigma, 2)) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S * normCDF(d1) - X * std::exp(-r * T) * normCDF(d2);
}

// Heston Model
double hestonModel(double S_0, double X, double T, double r, double kappa, double theta, double sigma_v, double rho, double v_0, int num_simulations = 10000, int num_steps = 100) {
    br::mt19937 gen;
    br::normal_distribution<> d(0, 1);

    double dt = T / num_steps;
    double option_payoff_sum = 0.0;

    for (int i = 0; i < num_simulations; ++i) {
        double S_t = S_0;
        double v_t = v_0;

        for (int j = 0; j < num_steps; ++j) {
            double z1 = d(gen);
            double z2 = rho * z1 + std::sqrt(1 - std::pow(rho, 2)) * d(gen);
            S_t += r * S_t * dt + std::sqrt(std::max(v_t, 0.0)) * S_t * z1 * std::sqrt(dt);
            v_t += kappa * (theta - v_t) * dt + sigma_v * std::sqrt(std::max(v_t, 0.0)) * z2 * std::sqrt(dt);
        }

        option_payoff_sum += std::max(S_t - X, 0.0);
    }

    return std::exp(-r * T) * (option_payoff_sum / num_simulations);
}

// Merton Jump Diffusion Model
double mertonJumpDiffusion(double S_0, double X, double T, double r, double sigma, double lambda_jump, double m_jump, double delta_jump, int num_simulations = 10000, int num_steps = 100) {
    br::mt19937 gen;
    br::normal_distribution<> normal_dist(0, 1);
    br::poisson_distribution<> poisson_dist(lambda_jump * T / num_steps);

    double dt = T / num_steps;
    std::vector<double> S_t(num_simulations, S_0);
    double option_payoff_sum = 0.0;

    for (int step = 0; step < num_steps; ++step) {
        for (int sim = 0; sim < num_simulations; ++sim) {
            double z = normal_dist(gen);
            int jumps = poisson_dist(gen);
            double jump_sum = 0.0;
            for (int j = 0; j < jumps; ++j) {
                jump_sum += m_jump + delta_jump * z;
            }
            S_t[sim] *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * z + jump_sum);
        }
    }

    for (double s : S_t) {
        option_payoff_sum += std::max(s - X, 0.0);
    }

    return std::exp(-r * T) * (option_payoff_sum / num_simulations);
}

int main() {
    // Parameters for the Black-Scholes model
    double S = 100;    // Current stock price
    double X = 110;    // Strike price
    double T = 1;      // Time to expiration in years
    double r = 0.03;   // Risk-free interest rate
    double sigma = 0.25;// Volatility

    // Black-Scholes Price
    double bsPrice = blackScholes(S, X, T, r, sigma);
    std::cout << "Black-Scholes Price: " << bsPrice << std::endl;

    // Parameters for the Heston model
    double kappa = 2.0;    // Mean reversion rate
    double theta = 0.06;   // Long-term mean of volatility
    double sigma_v = 0.2;  // Volatility of volatility
    double rho = -0.5;     // Correlation between the asset return and volatility
    double v_0 = 0.05;     // Initial variance

    // Heston Model Price
    double hestonPrice = hestonModel(S, X, T, r, kappa, theta, sigma_v, rho, v_0);
    std::cout << "Heston Model Price: " << hestonPrice << std::endl;

    // Parameters for the Merton Jump Diffusion model
    double lambda_jump = 0.75;  // Jump intensity
    double m_jump = -0.06;      // Mean jump size
    double delta_jump = 0.12;   // Jump-size volatility

    // Merton Jump Diffusion Price
    double mertonPrice = mertonJumpDiffusion(S, X, T, r, sigma, lambda_jump, m_jump, delta_jump);
    std::cout << "Merton Jump Diffusion Price: " << mertonPrice << std::endl;

    return 0;
}
