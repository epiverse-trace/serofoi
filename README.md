
## *serofoi* <img src="man/figures/serofoi-logo.png" align="right" width="120" />

An R package to estimates the *Force-of-Infection* of a given pathogen
from population based sero-prevalence studies on a Bayesian framework.

<!-- <!-- badges: start -->

–\>
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) -->
<!-- [![R-CMD-check](https://github.com/epiverse-trace/readepi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/epiverse-trace/readepi/actions/workflows/R-CMD-check.yaml) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/epiverse-trace/readepi/branch/main/graph/badge.svg)](https://app.codecov.io/gh/epiverse-trace/readepi?branch=main) -->
<!-- [![lifecycle-concept](https://raw.githubusercontent.com/reconverse/reconverse.github.io/master/images/badge-concept.svg)](https://www.reconverse.org/lifecycle.html#concept)  -->
<!-- <!-- badges: end --> –\>

## Installation

You can install the **development version** of `serofoi` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("TRACE-LAC/serofoi")
```

## Quick start

These examples illustrate some of the current functionalities:

# The package provides an example dataset of a serosurvey `mydata`

``` r

head(mydata)
#>           survey total counts age_min age_max year_init year_end tsur country
#> X.575 COL-035-18     2      0       1       1      2007     2007 2007     COL
#> X.576 COL-035-18     1      0       2       2      2007     2007 2007     COL
#> X.577 COL-035-18    13      2       4       4      2007     2007 2007     COL
#> X.578 COL-035-18    25      5       5       5      2007     2007 2007     COL
#> X.579 COL-035-18    17      0       6       6      2007     2007 2007     COL
#> X.580 COL-035-18    20      4       7       7      2007     2007 2007     COL
#>          test antibody
#> X.575 ELISA .      IgG
#> X.576 ELISA .      IgG
#> X.577 ELISA .      IgG
#> X.578 ELISA .      IgG
#> X.579 ELISA .      IgG
#> X.580 ELISA .      IgG

# Prepare the data for using serofoi

data_test <- prepare_data(mydata)


# Current version of the package runs three different models of the FoI

# Constant Force-of-Infection with a binomial distribution
model_0 <- run_model(model_data = data_test,
                     model_name = "constant_foi_bi",
                     n_iters = 500, 
                     n_thin = 2)
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] "serofoi model constant_foi_bi finished running ------"
#>            [,1]             
#> model_name "constant_foi_bi"
#> dataset    "COL-035-18"     
#> country    "COL"            
#> year       "2007"           
#> test       "ELISA ."        
#> antibody   "IgG"            
#> n_sample   "212"            
#> n_agec     "27"             
#> n_iter     "500"            
#> elpd       "-30.85"         
#> se         "4.9"            
#> converged  "Yes"

# Time-varying Force-of-Infection with a prior normal-binomial distribution
model_1 <- run_model(model_data = data_test,
                     model_name = "continuous_foi_normal_bi",
                     n_iters = 500, 
                     n_thin = 2)
#> Warning in readLines(file, warn = TRUE): incomplete final line found on
#> '/Library/Frameworks/R.framework/Versions/4.2/Resources/library/serofoi/extdata/stanmodels/continuous_foi_normal_bi.stan'
#> Warning in readLines(file, warn = TRUE): incomplete final line found on
#> '/Library/Frameworks/R.framework/Versions/4.2/Resources/library/serofoi/extdata/stanmodels/continuous_foi_normal_bi.stan'
#> Warning: There were 4 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.2, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
#> [1] "serofoi model continuous_foi_normal_bi finished running ------"
#>            [,1]                      
#> model_name "continuous_foi_normal_bi"
#> dataset    "COL-035-18"              
#> country    "COL"                     
#> year       "2007"                    
#> test       "ELISA ."                 
#> antibody   "IgG"                     
#> n_sample   "212"                     
#> n_agec     "27"                      
#> n_iter     "500"                     
#> elpd       "-29.74"                  
#> se         "5.25"                    
#> converged  NA

# Time-varying Force-of-Infection with a prior normal-log distribution
model_2 <- run_model(model_data = data_test,
                     model_name = "continuous_foi_normal_log",
                     n_iters = 500, 
                     n_thin = 2)
#> Warning in readLines(file, warn = TRUE): incomplete final line found on
#> '/Library/Frameworks/R.framework/Versions/4.2/Resources/library/serofoi/extdata/stanmodels/continuous_foi_normal_log.stan'
#> Warning in readLines(file, warn = TRUE): incomplete final line found on
#> '/Library/Frameworks/R.framework/Versions/4.2/Resources/library/serofoi/extdata/stanmodels/continuous_foi_normal_log.stan'
#> Warning: There were 104 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.33, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] "serofoi model continuous_foi_normal_log finished running ------"
#>            [,1]                       
#> model_name "continuous_foi_normal_log"
#> dataset    "COL-035-18"               
#> country    "COL"                      
#> year       "2007"                     
#> test       "ELISA ."                  
#> antibody   "IgG"                      
#> n_sample   "212"                      
#> n_agec     "27"                       
#> n_iter     "500"                      
#> elpd       "-29.7"                    
#> se         "5.45"                     
#> converged  "Yes"

# The following function can be used to visualize the sero-prevalence data with its corresponding binomial confidence interval before fitting to a model
plot_seroprev(data_test, size_text = 6)
```

<img src="man/figures/README-ex-1.png" width="100%" />

``` r

# For each model, you can use ploting functions such as


# You can compare these three models based on convergence, elpd and p-values
comp_table <- get_comparison_table(
  model_objects_list = c(m0 = model_0,
                         m1 = model_1,
                         m2 = model_2))
#> [1] "number of converged models = 2"
```

### Lifecycle

This package is currently a *concept*, as defined by the [RECON software
lifecycle](https://www.reconverse.org/lifecycle.html). This means that
essential features and mechanisms are still being developed, and the
package is not ready for use outside of the development team.

### Contributions

Contributions are welcome via [pull
requests](https://github.com/TRACE-LAC/serofoi/pulls).

Contributors to the project include:

- [Zulma M. Cucunubá](https://github.com/zmcucunuba) (author)

- \[Nicolás Tórres\] (author)

- \[Benjamin Lambert\] (author)

- \[Pierre Nouvellet\] (author)

- [Miguel Gamez](https://github.com/megamezl) (contributor)

- [Geraldine Gómez](https://github.com/megamezl) (contributor)

- [Jaime A. Pavlich-Mariscal](https://github.com/jpavlich) (contributor)

### Code of Conduct

Please note that the linelist project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
