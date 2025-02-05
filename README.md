
## *serofoi*: force-of-infection from population based serosurveys with age-disaggregated data <img src="man/figures/logo.png" align="right" width="130"/>

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/epiverse-trace/serofoi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/epiverse-trace/serofoi/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/epiverse-trace/serofoi/branch/dev/graph/badge.svg)](https://app.codecov.io/gh/epiverse-trace/serofoi/tree/dev/R?displayType=list)
[![lifecycle-maturing](https://raw.githubusercontent.com/reconverse/reconverse.github.io/master/images/badge-maturing.svg)](https://www.reconverse.org/lifecycle.html#maturing)

<!-- badges: end -->

*serofoi* is an R package to estimate the *Force-of-Infection (FoI)* of
a given pathogen from age-disaggregated population-based cross-sectional
serosurveys, using a Bayesian framework. The package provides a set of
features for assessing model fitting, convergence and visualisation.

*serofoi* relies on the
[`rstan`](https://mc-stan.org/users/interfaces/rstan) package, which
provides an R interface for the Stan programming language for
statistical Bayesian modelling. Particularly, *serofoi* relies on the
use of a *Hamiltonian Monte Carlo (HMC)* algorithm implemented by *Stan
for Markov chain Monte Carlo (MCMC)* sampling. The implemented methods
are outlined in ([Cucunubá et al. 2017](#ref-cucunubá2017)) and
([Carrera et al. 2020](#ref-carrera2020)) (see [FoI
Models](https://epiverse-trace.github.io/serofoi/articles/foi_models.html)
for further details). A compelling mathematical treatment of the
implemented serocatalytic models can be found in ([Kamau et al.
2025](#ref-kamau2025)).

*serofoi* is part of the [Epiverse
Initiative](https://data.org/initiatives/epiverse/).

## Installation

You can install the **development version** of *serofoi* from
[GitHub](https://github.com/epiverse-trace/serofoi) running:

``` r
if(!require("pak")) install.packages("pak")
pak::pak("epiverse-trace/serofoi")
```

or:

``` r
if(!require("remotes")) install.packages("remotes")
remotes::install_github("epiverse-trace/serofoi")
```

## Quick start

*serofoi* provides some minimal serosurvey datasets that can be used to
test out the package. For instance, the dataset `chagas2012` contains
seroprevalence measures of IgG antibodies against Trypanosoma cruzi
infection corresponding to a serological survey conducted in Colombia
during 2012 on a rural indigenous community that is known to present
long-term endemic transmission

``` r
# Load example dataset chagas2012 included with the package
data(chagas2012)
head(chagas2012, 5)
#>   survey_year n_sample n_seropositive age_min age_max
#> 1        2012       34              0       1       1
#> 2        2012       25              0       2       2
#> 3        2012       35              1       3       3
#> 4        2012       29              0       4       4
#> 5        2012       36              0       5       5
```

A visualisation of the serological data can be obtained using the
function `plot_serosurvey`:

``` r
plot_serosurvey(chagas2012, bin_serosurvey = TRUE, size_text = 15)
```

<img src="man/figures/README-data_test-1.png" width="50%" style="display: block; margin: auto;" />

A constant force-of-infection model can easily be implemented by means
of `fit_serodemol`:

``` r
seromodel <- fit_seromodel(serosurvey = chagas2012)
```

For further details on how to visualise the results and other available
models, please refer to the [online
documentation](https://epiverse-trace.github.io/serofoi/).

### Contributions

Contributors to the project include:

- [Zulma M. Cucunubá](https://github.com/zmcucunuba) (author,
  maintainer)

- [Nicolás Torres](https://github.com/ntorresd) (author)

- [Ben Lambert](https://ben-lambert.com/about/) (author)

- [Pierre Nouvellet](https://github.com/pnouvellet) (author)

- [Geraldine Gómez](https://github.com/megamezl) (contributor)

- [Jaime A. Pavlich-Mariscal](https://github.com/jpavlich) (contributor)

- [Miguel Gamez](https://github.com/megamezl) (contributor)

- [Hugo Gruson](https://github.com/Bisaloo) (contributor)

- [David Santiago Quevedo](https://github.com/davidsantiagoquevedo)
  (contributor)

- [Everlyn Kamau](https://github.com/ekamau) (contributor)

- [Richard Creswell](https://github.com/rccreswell) (contributor)

- [Sumali Bajaj](https://github.com/sumalibajaj) (contributor)

## Package vignettes

More details on how to use *serofoi* can be found in the [online
documentation as package
vignettes](https://epiverse-trace.github.io/serofoi/), under “Articles”.

## Help

To report a bug please open an
[issue](https://github.com/epiverse-trace/serofoi/issues/new/choose).

## Contribute

Contributions to *serofoi* are welcomed. Please follow the [package
contributing
guide](https://github.com/epiverse-trace/serofoi/blob/main/.github/CONTRIBUTING.md).

## Code of conduct

Please note that the *serofoi* project is released with a [Contributor
Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

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

<div id="ref-kamau2025" class="csl-entry">

Kamau, Everlyn, Junjie Chen, Sumali Bajaj, Nicolas Torres, Richard
Creswell, Jaime A Pavlich-Mariscal, Christl Donnelly, Zulma Cucunuba,
and Ben Lambert. 2025. “The Mathematics of Serocatalytic Models with
Applications to Public Health Data.” *medRxiv*, 2025–01.

</div>

</div>
