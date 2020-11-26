#' produce an MCMC-based EpiEstim R[t] estimate on the 29-day series
#' @import EpiEstim
#' @import shinytoastr
#' @param aout instance of 'Arima_sars2pack' S3 class
#' @note Simply uses the tsfull and dates29 components to call EpiEstim::estimate_R
#' @examples
#' nyd = nytimes_state_data()
#' ny = Arima_by_state(nyd, state.in="New York")
#' ee = est_Rt(ny)
#' od = options()$digits
#' options(digits=3)
#' tail(ee$R)
#' options(digits=od)
#' @export
est_Rt = function(aout) {
 #source(system.file("MCMCsupport/code.R", package="sars2death"))
 run_ee(aout)
}

#' basic support, formerly completely hidden in inst, uses Gamma-based MCMC for SI
#' @param arima_out instance of Arima_sars2pack
run_ee = function(arima_out, ...) {
 requireNamespace("EpiEstim")
 data(si_cens_G, package="sars2death")
 adf = data.frame(I=as.numeric(arima_out$tsfull), dates=arima_out$dates29)
 EpiEstim::estimate_R(adf, method="si_from_sample", si_sample=si_cens_G$si_sample,
     config=EpiEstim:::make_config(list(n2=50, seed=1)))
}

