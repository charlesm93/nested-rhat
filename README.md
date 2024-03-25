# Nested $\widehat R$: supplemental code

Supplemental code for the [preprint](https://arxiv.org/abs/2110.13017) _Nested_ $\widehat R$: _Assessing the convergence of Markov chains Monte Carlo when running many short chains_ by Charles Margossian, Matt Hoffman, Pavel Sountsov, Lionel Riou-Durand, Aki Vehtari, and Andrew Gelman.

The heavy-duty experiments are done in the `jupyter` notebooks. We run `TensorFlow Probability` built on `jax`. The code runs on both CPU and GPU, however we strongly recommend using GPU when running a large number of Markov chains.

## Figures and how to reproduce them

* Figure 1: Rosenbrock distribution example. Generated by `rosenbrock_example.ipynb`.
* Figure 3: Variance of Monte Carlo estimators with constrained and naive superchains. Generated by `nonstationary_variance.ipynb`.
* Figure 4: Variance of $\widehat B_{\mathfrak n} / \widehat W_{\mathfrak n}$. Generated by `variance_simulation.R`
* Figures 5-10: Numerical experiments. The data for the figures are generated by `nRhat-convergence.ipynb`, with `pk_simulation.ipynb` as a supporting file. The data is then read by `analyze_results.ipynb`, with which we generate the figures.
* Figure 11 (Appendix): Threshold for $(\delta, \delta')$-reliability. Generated by `langevin_simulation.R`.  



## Details for `nRhat_convergence.ipynb`

At the top of the notebook, you can specify the control variables `num_chains` and `num_super_chains` (the number of distrinct initializations), as well as whether or not to use naive superchains (`naive_super_chains = True`). The notebook applies these control parameters to all 6 models in the Numerical Experiment section and saves the data. To generate the plots in the paper, we use `num_chains=2048` and vary `num_super_chains` (2, 8, 16, 64, 256, 1024). Each run takes 1 - 2 hours on a GPU.
