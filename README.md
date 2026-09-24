# Stochastic Processes Lab

Numerical experiments and visualizations for stochastic processes, probability, and Monte Carlo methods.

This undergraduate mathematics project will connect **mathematical theory → numerical simulation → empirical verification**. It currently contains only the project structure; no experiments or results have been completed.

## Planned experiments

1. [Random Walk](notebooks/01_random_walk.ipynb)
2. [Poisson Process](notebooks/02_poisson_process.ipynb)
3. [Markov Chain](notebooks/03_markov_chain.ipynb)
4. [Ornstein–Uhlenbeck Process](notebooks/04_ornstein_uhlenbeck.ipynb)
5. [Monte Carlo Methods](notebooks/05_monte_carlo.ipynb)

Each notebook is reserved with the same eight sections: Problem, Mathematical Background, Theoretical Result, Numerical Experiment, Visualization, Comparison with Theory, Interpretation, and Limitations. Future experiments should state the mathematical question, derive or cite the theoretical result, run a numerical simulation with a fixed random seed, and explain how the results compare with theory.

## Project layout

- `notebooks/`: explanations and visualizations
- `src/`: reusable simulation and analysis functions
- `tests/`: tests for important functions
- `figures/`: generated figures, once experiments exist
- `notes/`: mathematical notes and derivations

Dependencies planned for later experiments are listed in `requirements.txt`. No numerical results are claimed at this stage.
