# serofoi (Development version)

* Updated the stan code to use new syntax.

# serofoi 1.0.2

## New features

* **Serological surveys simulation functionalities**
  * Introduce enhanced data simulation functionalities to support a broader range of use cases for serosurvey analysis.
  * Add a dedicated vignette for data simulation to guide users on simulating data for serosurveys.

* **Enables the implementation of a wide variety of new models**
  * Add forward random walk age-varying models with uniform and normal priors.
  * Improves previous implementation of forward random walk time-varying models.
  * Enables prior distributions specifications by means of package specific functions like `sf_normal()` and `sf_uniform()`.
  * Support for seroreversion rate estimation for all models with uniform and normal priors.

## Breaking changes

* Replaced old modelling and visualization functions, making Bayesian model selection and specification more flexible.
* New models and functionalities include constant, time-varying, and age-varying models, as well as options for estimating seroreversion rates.
* Updated the R-hat convergence threshold for model convergence diagnostics to R-hat<1.01,
following [Vehtari et al. (2021)](https://projecteuclid.org/journals/bayesian-analysis/volume-16/issue-2/Rank-Normalization-Folding-and-Localization--An-Improved-R%cb%86-for/10.1214/20-BA1221.full).

## Minor changes

* Remove the `simdata_*` datasets from the package and replaced them with code-based simulation in vignettes.
* Add missing examples to exported functions
* Add missing documentation entries

## Internal changes

* **Unit testing**
  * Added unit tests for key package functionalities, including:
    - `fit_seromodel()` function (#213).
    - Visualization functionalities (#215).
    - `build_stan_data` and related functions (#232).
    - Validation functions (#235).
  * Test coverage increased to 100%

* **Refactorization of simulation examples**
  * Replaced static simulated datasets in tests and vignettes with dynamic examples using data simulation functions.

* **Documentation improvements**
  * Enhanced documentation for new functionalities and updated vignettes to reflect recent changes.

# serofoi 0.1.0

## New features

* Data simulation functions from time-varying or age-varying _force of infection_ trends. The following is an example to simulate from a constant (in time) _force of infection_:
```
foi_sim_constant <- rep(0.02, 50)

serodata_constant <- generate_sim_data(  
  sim_data = data.frame(  
    age = seq(1, 50),  
    tsur = 2050),  
  foi = foi_sim_constant,  
  sample_size_by_age = 5  
)  
```
To generate grouped serosurveys the function `group_sim_data` can be used:
```
serodata_constant <- group_sim_data(serodata_constant, step = 5)
```

## Breaking changes

* Simplifies `fit_seromodel` output

  * Before, the output of `fit_seromodel` was a list:
```
    seromodel_object <- list(
      fit = fit,
      seromodel_fit = seromodel_fit,
      serodata = serodata,
      serodata = serodata,
      stan_data = stan_data,
      ...
    )
```
  * Now, the output is a `stan_fit` object as obtained from [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html). Because of this, **some plotting functionalities now require `serodata` as an input**.

* Initial prior distribution parameters `foi_location` and `foi_scale` can be specified explicitly in `fit_seromodel`:
```
seromodel <- fit_seromodel(
  serodata,
  foi_model = "tv_normal",
  foi_location = 0,
  foi_scale = 1
)
```
Depending on the selected model `foi_model`, the meaning of the parameters change. For the `tv_normal_log` model these parameters must be in logarithmic scale; the recommended usage is:
```
seromodel <- fit_seromodel(
  serodata,
  foi_model = "tv_normal_log",
  foi_location = -6,
  foi_scale = 4
)
```

* **Chunks structure specification is now possible**
  * Before, the models estimated one value of the _force of infection_ per year in the time spanned by the serosurvey:
```
data(chagas2012)
serodata <- prepare_serodata(chagas2012)
seromodel <- fit_seromodel(serodata, foi_model = "tv_normal")
```
![image](https://github.com/epiverse-trace/serofoi/assets/45337127/3ab8e761-d92b-4d10-8897-a8c6d6add854) 

  * Now, the amount of _force of infection_ values estimated by the models depend on the specified chunk structure. This can either be specified by size:
```
seromodel <- fit_seromodel(serodata, foi_model = "tv_normal", chunk_size = 10)
```
![image](https://github.com/epiverse-trace/serofoi/assets/45337127/b70e2315-64b5-4cbb-b770-85b6f27175e8)

or explicitly:
```
chunks <- rep(c(1, 2, 3, 4, 5), c(10, 10, 15, 15, max(serodata$age_mean_f)-50))
seromodel <- fit_seromodel(serodata, foi_model = "tv_normal", chunks = chunks)
```
![image](https://github.com/epiverse-trace/serofoi/assets/45337127/2cb998db-b86b-4c2d-9693-4a683d3a1267)

* **Deprecate `run_seromodel`**. Initially this function was intended to be a handler for `fit_seromodel` for cases when the user may need to implement the same model to multiple independent serosurveys; now we plan to showcase examples of this using the current functionalities of the package (to be added in future versions to the vignettes). 

## Minor changes

* **Refactorization of the visualization module**
  1. `plot_seroprev` allows for data binning (age group manipulation) by means of parameters `bin_data=TRUE` and `bin_step`.
  2. Automatic selection of `ymin` and `ymax` aesthetics plotting functions (with the exception of `plot_rhats`).
  3. Correct input validation

* Remove duplicated data in `veev2012` dataset

## Internal changes

* Remove large files from git history (see #77).

* Added input validation for the following functions:
  * `prepare_serodata`
  * `generate_sim_data`
  * `get_age_group`
  * `fit_seromodel`
  * `extract_seromodel_summary`
  * `plot_seroprev`
  * `plot_seroprev_fitted`
  * `plot_foi`
  * `plot_seromodel`

* Unit testing:
  * Separate modelling tests by model
  * Use of  `dplyr::near` to test models statistical validity
  * Add tests for data simulation functions

* Update package template in accordance to [{packagetemplate}](https://github.com/epiverse-trace/packagetemplate)

# serofoi 0.0.9

This release of _**serofoi**_, includes the following:

1. Implementation of package modules: Incorporates data preparation, modelling, and visualization modules, they enable efficient handling of data, perform statistical modelling, and generate visual representations of the results.
2. Documentation: It consists of vignettes, a website, and uses cases that provide detailed instructions on how to use the package effectively.
3. Implementation of 3 models for calculating the Force-of-Infection (FoI): The first model is the constant or endemic model, which assumes a stable FoI over time. The second and third models are time-varying, with the normal FoI model representing a slow change in FoI and the normal-log FoI model representing a fast epidemic change in FoI.
4. Definition of coverage test to assurance the quality of the package.

Overall, this release introduces essential package functionality, comprehensive documentation, various FoI models, and a coverage test, enabling users to analyse seroprevalence data and calculate the Force-of-infection.
