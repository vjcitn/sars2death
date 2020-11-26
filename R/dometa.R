#' obtain input to rmeta meta.summaries, metaplot
#' @import rmeta
#' @import parallel
#' @param nyd tibble, data source assumed to be nytimes_state_data()
#' @param opt_parms list, such as that produced by `min_bic_all_states()`
#' @param \dots passed to `Arima_by_state`
#' @examples
#' nyd = nytimes_state_data()
#' data(min_bic_2020_11_26_death)
#' m1 = run_meta(nyd, opt_parms=min_bic_2020_11_26_death) # must be relatively current
#' rmeta::meta.summaries(m1$drifts, m1$se.drifts)
#' names(m1$drifts) = gsub(".drift", "", names(m1$drifts))
#' nyind = which(names(m1$drifts) %in% c("New York", "New Jersey"))
#' rmeta::meta.summaries(m1$drifts[-nyind], m1$se.drifts[-nyind])
#' o = order(m1$drifts)
#' rmeta::metaplot(m1$drifts[o], m1$se.drifts[o], labels=names(m1$drifts)[o], cex=.7, 
#'   xlab="Infection velocity (CHANGE in number of confirmed cases/day)", ylab="State")
#' segments(rep(-350,46), seq(-49,-4), rep(-50,46), seq(-49,-4), lty=3, col="gray")
#' @export
run_meta = function(nyd, opt_parms, ...) {
 allst = contig_states_dc()
 nc = parallel::detectCores()-1
 options(mc.cores=nc)
 allarima = lapply(allst, function(x) {
   parms = opt_parms[[x]]
   Arima_by_state(nyd, x, ARorder=parms$opt["ARord"], MAorder=parms$opt["MAord"], ...)
   })
 names(allarima) = allst
 drifts = sapply(allarima, function(x) coef(x$fit)["drift"])
 searima = function(a) sqrt(a$fit$var.coef["drift", "drift"])
 se.drifts = sapply(allarima, searima)
 list(drifts=drifts, se.drifts=se.drifts)
}

#' county within state metaanalysis of incidence velocity
#' @param nytc nytimes_county_data() output
#' @param state.in character(1)
#' @param verbose logical(1) will print '.' after each county's optimal AR/MA orders found
#' @param \dots passed to `Arima_by_county`
#' @note The search for optimal AR/MA orders may throw warnings.
#' @examples
#' nytc = nytimes_county_data()
#' m1 = run_meta_county(nytc)
#' names(m1$drifts) = gsub(".drift", "", names(m1$drifts))
#' rmeta::meta.summaries(m1$drifts, m1$se.drifts)
#' o = order(m1$drifts)
#' rmeta::metaplot(m1$drifts[o], m1$se.drifts[o], labels=names(m1$drifts)[o], cex=.7, 
#'   xlab="Infection velocity (CHANGE in number of confirmed cases/day)", ylab="County")
#' @export
run_meta_county = function(nytc, state.in="Massachusetts", verbose=TRUE, ...) {
 allc = unique(nytc$county[which(nytc$state==state.in)])
 nc = parallel::detectCores()-1
 options(mc.cores=nc)
 allarima = lapply(allc, function(x) {
   if (verbose) cat(".")
   Arima_by_county(nytc, state.in=state.in, county.in=x, ARorder=NULL, MAorder=NULL, ...)
   })
 ng = sapply(allarima, inherits, "try-error")
 dropped = NULL
 names(allarima) = allc
 if (sum(ng)>0) {
    message("could not find BIC-optimal ARMA orders for some counties")
    print(dropped <- allc[which(ng)])
    message("dropping data on these")
    allarima = allarima[-which(ng)]
    }
 drifts = sapply(allarima, function(x) coef(x$fit)["drift"])
 searima = function(a) sqrt(a$fit$var.coef["drift", "drift"])
 se.drifts = sapply(allarima, searima)
 ans = list(drifts=drifts, se.drifts=se.drifts, dropped=dropped)
 class(ans) = "meta_county_sars2death"
 ans
}

#' @export
plot.meta_county_sars2death = function(x, y, ..., cex.in=.7) {
 names(x$drifts) = gsub(".drift", "", names(x$drifts))
 o = order(x$drifts)
 rmeta::metaplot(x$drifts[o], x$se.drifts[o], labels=names(x$drifts)[o], cex=cex.in,
   xlab="Infection velocity (CHANGE in number of confirmed cases/day)", ylab="County")
}
