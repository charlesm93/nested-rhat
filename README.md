# Nested $\widehat R$: supplemental code

Supplemental code for the [preprint](https://arxiv.org/abs/2110.13017) Nested $\widehat R$: Assessing the convergence of Markov chains Monte Carlo when running many short chains by Charles Margossian, Matt Hoffman, Pavel Sountsov, Lionel Riou-Durand, Aki Vehtari, and Andrew Gelman.

The Python Colab notebooks can be used to reproduce our numerical experiments. Most of these notebooks are forked from [google-research/nested-rhat](https://github.com/google-research/google-research/tree/master/nested_rhat), which contains the code for the first version of our preprint.

All Python experiments are run using Colab Pro. We recommend using a GPU and the high ram setting.

Some of our less computationally intensive experiments are written in R.

## Figures and how to reproduce them

* Figure 2: Example on Banana problem. Use `rhat_locker.ipynb`.
* Figure 3: Variance of Monte Carlo estimators. Use `Nested_R_hat.ipynb`.
* Figures 4, 11, 12: Variance of $\mathfrak n \widehat B / \mathfrak n \widehat W$. Use `variance_simulation.R`.
* Figures 5, 6, 7, 8, 13, 14: Numerical experiments (Section 5). Use `nRhat_convergence.ipynb`, see details.
* Figure 9. $(\delta, \delta')$-reliability. Use `Bias-Variance-MCMC.ipynb`.

## Details for `nRhat_convergence.ipynb`

The more computationally intensive experiments happen in `nRhat_convergence.ipynb`. Unfortunately, all figures cannot be generated in one go without exhausting the ram on a GPU. It is therefore recommended to run the experiments, save the results, and then generate the figures. The notebook saves experimental results to your Google Drive, see last cell before _Application to models_, 
```
# Packages to save files on the google drive
from google.colab import files
from google.colab import drive

drive.mount('/content/gdrive')

dir_name = "gdrive/MyDrive/nRhat_experiments/"
```
which can be adjusted to your setting.

After running the cells under `Application to models/Setup`, the user can choose which example to experiment with. We'll here consider the _Banana_ problem. First run the cells under the _Banana_ header. The experiments under the header _$\mathfrak n \widehat R$ diagnostic_ were not included in the paper, and you may skip ahead to _Adaptive warmup length_. MCMC is run on the Banana target with recordings of $\mathfrak n \widehat R$ and the squared error recorded at the end of each warmup window (see `window_array`). Be sure to check the setting of MCMC (number of superchains, samples, etc). The experiment is repeated `num_seeds` times. The results are then saved on the mounted drive.

Under the _Model synthesis_ header, experimental results can be read to generate plots.
