#' Extract specified loo estimate
#'
#' @inheritParams extract_central_estimates
#' @param par_loo_estimate Name of the loo estimate to be extracted.
#' @param loo_estimate_digits Number of loo estimate digits
#' @return Text summarising specified loo estimate
#' @export
summarise_loo_estimate <- function(
    seromodel,
    par_loo_estimate = "elpd_loo",
    loo_estimate_digits = 2
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  loo_fit <- loo::loo(
    seromodel,
    pars = c(parameter_name = "log_likelihood")
  )
  loo_estimate <- loo_fit$estimates[par_loo_estimate, ] %>%
    round(loo_estimate_digits)

  loo_estimate_summary <- paste0(loo_estimate[1], "(se=", loo_estimate[2], ")")

  return(loo_estimate_summary)
}

#' Summarise central estimate
#'
#' @inheritParams extract_central_estimates
#' @param central_estimate_digits Number of central estimate digits
#' @return Text summarising specified central estimate
#' @export
summarise_central_estimate <- function(
    seromodel,
    serosurvey,
    alpha,
    par_name = "seroreversion_rate",
    central_estimate_digits = 2
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  central_estimates <- signif(
    extract_central_estimates(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = par_name
    ),
    digits = 2
  )

  central_estimate_summary <- paste0(
    central_estimates$median,
    "(", 100 * (1 - alpha), "% CI, ",
    central_estimates$lower, "-",
    central_estimates$upper, ")"
  )

  return(central_estimate_summary)
}

#' Summarise specified model
#'
#' @inheritParams extract_central_estimates
#' @inheritParams summarise_loo_estimate
#' @inheritParams summarise_central_estimate
#' @return A list summarising the specified model
#' \describe{
#'  \item{`model_name`}{Name of the model}
#'  \item{`elpd`}{elpd and its standard deviation}
#'  \item{`foi`}{Estimated foi with credible interval (for 'constant' model)}
#'  \item{`foi_rhat`}{foi rhat value (for 'constant' model)}
#'  \item{`seroreversion_rate`}{Estimated seroreversion rate}
#'  \item{`seroreversion_rate_rhat`}{Seroreversion rate rhat value}
#' }
#' @export
summarise_seromodel <- function(
    seromodel,
    serosurvey,
    alpha = 0.05,
    par_loo_estimate = "elpd_loo",
    loo_estimate_digits = 1,
    central_estimate_digits = 2,
    rhat_digits = 2
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  model_name <- seromodel@model_name
  summary_list <- list(model_name = model_name)

  loo_estimate_summary <- summarise_loo_estimate(
    seromodel = seromodel,
    par_loo_estimate = par_loo_estimate,
    loo_estimate_digits = loo_estimate_digits
  )

  summary_list[par_loo_estimate] <- loo_estimate_summary

  check_convergence <- NULL
  if (startsWith(model_name, "constant")) {
    foi_summary <- summarise_central_estimate(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = "foi",
      central_estimate_digits = central_estimate_digits
    )

    foi_rhat <- bayesplot::rhat(seromodel, "foi") %>%
      signif(rhat_digits)

    check_convergence <- append(
      check_convergence,
      foi_rhat < 1.01
    )

    summary_list <- append(
      summary_list,
      list(
        foi = foi_summary,
        foi_rhat = foi_rhat
      )
    )
  } else {
    rhats <- bayesplot::rhat(seromodel, "foi_vector")

    check_convergence <- append(
      check_convergence,
      all(rhats < 1.01)
    )
  }

  if (!endsWith(model_name, "no_seroreversion")) {
    seroreversion_rate_summary <- summarise_central_estimate(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = "seroreversion_rate",
      central_estimate_digits = central_estimate_digits
    )

    seroreversion_rate_rhat <- bayesplot::rhat(
      seromodel,
      "seroreversion_rate"
      ) %>%
      signif(rhat_digits)

    check_convergence <- append(
      check_convergence,
      seroreversion_rate_rhat < 1.01
    )

    summary_list <- append(
      summary_list,
      list(
        seroreversion_rate = seroreversion_rate_summary,
        seroreversion_rate_rhat = seroreversion_rate_rhat
      )
    )
  }

  if (all(check_convergence)) {
    summary_list["converged"] <- "yes"
  } else {
    summary_list["converged"] <- "no"
  }

  return(summary_list)
}
