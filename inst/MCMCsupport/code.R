
run_ee = function(arima_out, ...) {
 requireNamespace("EpiEstim")
 load(system.file("MCMCsupport/si_cens_G.rda", package="sars2app"))
 adf = data.frame(I=as.numeric(arima_out$tsfull), dates=arima_out$dates29)
 EpiEstim::estimate_R(adf, method="si_from_sample", si_sample=si_cens_G$si_sample, 
     config=EpiEstim:::make_config(list(n2=50, seed=1)))
}
