
## *serofoi*: force-of-infection from population based serosurveys with age-disagregated data <img src="man/figures/serofoi-logo.png" align="right" width="130"/>

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/TRACE-LAC/serofoi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TRACE-LAC/serofoi/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/TRACE-LAC/serofoi/branch/bugfixes-jaime/graph/badge.svg)](https://github.com/TRACE-LAC/serofoi?branch=bugfixes-jaime)
[![lifecycle-concept](https://raw.githubusercontent.com/reconverse/reconverse.github.io/master/images/badge-concept.svg)](https://www.reconverse.org/lifecycle.html#concept)

<!-- badges: end -->

***serofoi*** is an R package to estimates the *Force-of-Infection* of a
given pathogen from age-disaggregated population based sero-prevalence
studies, using a Bayesian framework.

***serofoi*** implements methods outlined in ([Cucunubá et al.
2017](#ref-cucunubá2017)) and ([Carrera et al. 2020](#ref-carrera2020))

***serofoi*** relies on
[`rstan`](https://mc-stan.org/users/interfaces/rstan)

***serofoi*** is part of the [Epiverse
Initiative](https://data.org/initiatives/epiverse/).

## Installation

You can install the **development version** of `serofoi` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
# remotes::install_github("TRACE-LAC/serofoi")
library(serofoi)
```

## Quick start

The package provides an example dataset of the observed serosurvey data,
`mydata`. This example is the basic entry for the package.

``` r
head(mydata)
#>       survey total counts age_min age_max year_init year_end tsur country
#> 1 COL-035-18     2      0       1       1      2007     2007 2007     COL
#> 2 COL-035-18     1      0       2       2      2007     2007 2007     COL
#> 3 COL-035-18    13      2       4       4      2007     2007 2007     COL
#> 4 COL-035-18    25      5       5       5      2007     2007 2007     COL
#> 5 COL-035-18    17      0       6       6      2007     2007 2007     COL
#> 6 COL-035-18    20      4       7       7      2007     2007 2007     COL
#>      test antibody
#> 1 ELISA .      IgG
#> 2 ELISA .      IgG
#> 3 ELISA .      IgG
#> 4 ELISA .      IgG
#> 5 ELISA .      IgG
#> 6 ELISA .      IgG
```

The function `prepare_data` will prepare the entry data for entering the
modelling functions. The seroprevalence *prepared data* can be
visualised with the `plot_seroprev` function. This function also plots
the binomial confidence interval of the observed data.

``` r

data_test <- prepare_data(mydata)

plot_seroprev(data_test, size_text = 15)
```

<img src="man/figures/README-data_test-1.png" width="100%" />

#### Current version of the package runs ***three*** different FoI models

The `run_model` function allows specifying the Bayesian model from *R*,
while running in the back from `rstan`. The number of iterations,
thinning, and other parameters can be customised.

<div class="alert alert-primary">

NOTE: Running the *serofoi* models for the first time on your local
computer make take a few minutes while the *rstan* code is compiled
locally. Afterwards, no further compilation is needed.

</div>

#### Model 1. Constant Force-of-Infection (endemic model)

For the *endemic model* a small number of iterations is enough for
achieving convergence, as it only fits one parameter (the constant FoI)
from a binomial distribution.

``` r
model_1 <- run_model(model_data = data_test,
                     model_name = "constant_foi_bi",
                     n_iters = 500, 
                     n_thin = 2)
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
```

#### Model 2. Time-varying Force-of-Infection (epidemic model)

For the *epidemic model,* a larger number of iterations is required for
achieving convergence, as it fits yearly FoI values from a binomial
distribution. The number of iterations required may depend on the number
of years, reflected by the difference between year of the serosurvey and
the maximum age-class sampled.

``` r
model_2 <- run_model(model_data = data_test,
                     model_name = "continuous_foi_normal_bi",
                     n_iters = 1500, 
                     n_thin = 2)
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
#> n_iter     "1500"                    
#> elpd       "-29.82"                  
#> se         "5.28"                    
#> converged  "Yes"
```

#### Model 3. Time-varying Force-of-Infection (fast epidemic model)

For the *fast* *epidemic model,* a larger number of iterations is
required for achieving convergence, compared to the previous models.

``` r
model_3 <- run_model(model_data = data_test,
                     model_name = "continuous_foi_normal_log",
                     n_iters = 1500, 
                     n_thin = 2)
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
#> n_iter     "1500"                     
#> elpd       "-29.73"                   
#> se         "5.53"                     
#> converged  "Yes"
```

For each model, the plot_model function generate a vertical arrange of
plots summarising the results of the model implementation plotting
functions. Crucially, it shows the (expected) log-predictive density
`elpd`, standard error `se`, and allows to check convergence based on
`R-hat` convergence diagnostics.

Also, the `plot_models_list` allows a visual a comparison of the models
based on the (expected) log-predictive density `elpd`, standard error
`se`, and allows to check convergence based on `R-hat` convergence
diagnostics.

For more detailed information and examples, please check the [online
documentation](https://epiverse-trace.github.io/serofoi/articles) as
package vignettes under Get Started.

### Contributions

Contributors to the project include:

- [Zulma M. Cucunubá](https://github.com/zmcucunuba) (author,
  maintainer)

- [Nicolás Tórres](https://github.com/ntorresd) (author)

- [Benjamin Lambert](https://ben-lambert.com/about/) (author)

- [Pierre Nouvellet](https://github.com/pnouvellet) (author)

- [Miguel Gamez](https://github.com/megamezl) (contributor)

- [Geraldine Gómez](https://github.com/megamezl) (contributor)

- [Jaime A. Pavlich-Mariscal](https://github.com/jpavlich) (contributor)

## Package vignettes

More details on how to use ***serofoi*** can be found in the [online
documentation as package
vignettes](https://epiverse-trace.github.io/serofoi/), under “Get
Started”.

## Help

To report a bug please open an
[issue](https://github.com/epiverse-trace/serofoi/issues/new/choose).

## Contribute

Contributions to ***serofoi*** are welcomed. Please follow the [package
contributing
guide](https://github.com/epiverse-trace/serofoi/blob/main/.github/CONTRIBUTING.md).

## Code of conduct

Please note that the ***serofoi*** project is released with a
[Contributor Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-carrera2020" class="csl-entry">

Carrera, Jean-Paul, Zulma M. Cucunubá, Karen Neira, Ben Lambert, Yaneth
Pittí, Jesus Liscano, Jorge L. Garzón, et al. 2020. “Endemic and
Epidemic Human Alphavirus Infections in Eastern Panama: An Analysis of
Population-Based Cross-Sectional Surveys.” *The American Journal of
Tropical Medicine and Hygiene* 103 (6): 2429–37.
<https://doi.org/10.4269/ajtmh.20-0408>.

</div>

<div id="ref-cucunubá2017" class="csl-entry">

Cucunubá, Zulma M, Pierre Nouvellet, Lesong Conteh, Mauricio Javier
Vera, Victor Manuel Angulo, Juan Carlos Dib, Gabriel Jaime Parra -Henao,
and María Gloria Basáñez. 2017. “Modelling Historical Changes in the
Force-of-Infection of Chagas Disease to Inform Control and Elimination
Programmes: Application in Colombia.” *BMJ Global Health* 2 (3):
e000345. <https://doi.org/10.1136/bmjgh-2017-000345>.

</div>

</div>
