# Option Pricing Models in C++

## Project Overview
This repository represents an experimental endeavor aimed at translating my previously completed financial analysis tools from Python to C++. The primary objective is to re-implement option pricing models, encompassing the Black-Scholes Model, the Heston Model, and the Merton Jump Diffusion Model. This project serves as a personal exploration, allowing me to assess the intricacies and performance disparities between C++ and Python within the realm of computational finance.

## Features
- **Black-Scholes Model**: Implementation of the Black-Scholes formula for option pricing in C++.
- **Heston Model**: Implementation of the Heston stochastic volatility model for option pricing, incorporating mean reversion and volatility of volatility.
- **Merton Jump Diffusion Model**: Implementation of the Merton Jump Diffusion model for option pricing, considering abrupt price changes with jumps.
- **Monte Carlo Simulation**: Utilizes Monte Carlo simulation to calculate option prices for more accurate financial modeling.

## Installation
To set up this C++ project, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/Gouldh/Option-Pricing-CPP.git
   ```
2. Navigate to the project directory:
   ```bash
   cd Option-Pricing-CPP
   ```
## Usage
To use the C++ option pricing models, you need a C++ development environment. You can compile and run the project using a C++ compiler such as `g++`. Here's how to compile and execute the code:

```bash
g++ -o option_pricing main.cpp -lboost_random
./option_pricing
```
Ensure that you have the `Boost` C++ Libraries installed, as the code uses `boost::random` for random number generation.

## License
This project is open-sourced under the MIT License. For more information, please refer to the `LICENSE` file.

**Author**: Hunter Gould
**Date**: 11/26/2023
