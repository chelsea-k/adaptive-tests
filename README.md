This repository contains scripts for reproducing the results from
Krantsevich et al., "Bayesian Decision Theory for Tree-Based
Adaptive Screening Tests with an Application to Youth Delinquency."

First, install the two packages in `adaptive-tests/packages`, which
are `bfa.mod` and `rpartMaxVPP`. `bfa.mod` is a modified 
version of the `bfa` package; the modifications allow for making
posterior draws indexed by draw number, and sampling from the 
conditional predictive distribution. `rpartMaxVPP` is a modified
version of the `rpart` package; the modifications allow for 
growing an `rpart` tree using the `maxvpp` criterion; for more
information see Krantsevich et al. 

To install these packages, after cloning the repo, in terminal enter 
`R CMD install <path_to_package_file>`. Or, in RStudio you can do 
`install.packages("<path_to_package_file>", type="source", repos=NULL)`
`<path_to_package_file>` is e.g. `packages/bfa.mod_0.4.tar.gz`.

Once the packages are installed, the scripts are generally supposed 
to be run in this order:

1) `simulate_data.R`: creates simulated item response data 

2) `draw_synth_XB.R`: creates synthetic data from the simulated data using
the Gaussian copula factor and XBART models.
`draw_synth_RF.R`: creates synthetic data from the simulated data using 
local perturbations and a Random Forest model.

3) `fit_CART_maxIPP.R`: fits a regression tree using the maxIPP and
pruning strategy described in the paper.
`fit_CART_maxDepth.R`: fits a regression and classification tree
using a maximum depth criterion.

4) `plot_Fig*.R` and `create_Table*.R`: reproduce the figures and 
tables from the paper using data and results generated in steps 1) 
through 3). Different figures and tables require data fit on all data
or only on training data; or fit to a subpopulation or to the whole
population. These settings are achieved through the boolean variables 
`out.of.sample` and `subpopulation` at the beginning of the `draw_synth_*.R`
and `fit_CART_*.R` scripts. To reproduce all plots and tables, run steps 1) 
through 3) with all four combinations (e.g. `out.of.sample = TRUE, FALSE`
and `subpopulation = TRUE, FALSE`), then run the plotting and table
scripts.

Note: `util_functions.R` contains several functions that are reused throughout
the plotting and table scripts.











