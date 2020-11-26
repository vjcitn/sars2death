#' generalize type of event modeled
#' @keywords internal
#event_tag = function() "confirmed"  # change to "deaths" if you want
event_tag = function() "deaths"  # change to "deaths" if you want


#' a two-dimensional grid of integers
#' @param n maximum value, will expand.grid(0:n,0:n) with colnames arorder, maorder
#' @export
grid_2d = function(n=5) expand.grid(arorder=seq(0,n), maorder=seq(0,n))

#' minimize BIC for ARIMA models with default differencing order 1 over choices of MA and AR order
#' @param src tibble with cumulative incidence like output of nytimes_state_data()
#' @param fullusa logical(1) if TRUE, use Arima_nation (src should be enriched_jhu_data())
#' @param state.in character(1) state name
#' @param county.in character(1) county name, defaults to NULL, but if supplied, focus on county-level data
#' @param parms two-column data frame with proposed AR order in column 1 and MA order in column 2
#' @param max_date character(1) or lubridate date from which to look back
#' @param lookback_days numeric(1)
#' @param \dots passed to `Arima_by_state` or `Arima_by_county`
#' @examples
#' nyd = nytimes_state_data()
#' mb = min_bic(nyd, state.in="New York")
#' names(mb)
#' nytc = nytimes_county_data()
#' mb2 = min_bic(nytc, state.in="Massachusetts", county.in="Norfolk")
#' names(mb2)
#' @export
min_bic = function(src, fullusa=FALSE, state.in="New York", county.in=NULL, parms=grid_2d(5), max_date=NULL, 
   lookback_days=29, ...) {
 nr = nrow(parms)
 if (!fullusa) bics = sapply(seq_len(nr), function(r) {
    if (is.null(county.in)) {
         ari = try(Arima_by_state(src, state.in=state.in, ARorder=parms$arorder[r],
          MAorder=parms$maorder[r], max_date=max_date, lookback_days=lookback_days, ...), silent=TRUE)
         if (!inherits(ari, "try-error")) return(ari$fit$bic)
         ari
         } else {
         ari = try(Arima_by_county(src, state.in=state.in, county.in=county.in, ARorder=parms$arorder[r],
          MAorder=parms$maorder[r], max_date=max_date, lookback_days=lookback_days, ...), silent=TRUE)
         if (!inherits(ari, "try-error")) return(ari$fit$bic)
         ari
         }
      })
 else {
       state.in = NULL
       bics = sapply(seq_len(nr), function(r) {
          ari = try(Arima_nation(src, ARorder=parms$arorder[r], MAorder=parms$maorder[r], 
               max_date=max_date, lookback_days=lookback_days, ...), silent=TRUE)
          if (!inherits(ari, "try-error")) return(ari$fit$bic)
          ari
      })
     }
 errs = sapply(bics, function(x) (inherits(x, "try-error") | is.na(x)))
 dr = which(errs)
 if (length(dr)>0) {
    bics=bics[-dr]
    parms = parms[-dr,]
    }
 ind = which.min(bics)
 md = max_date
 if (is.null(md)) md = max(src$date)
 ans = list(opt=c(ARord=parms$arorder[ind], MAord=parms$maorder[ind]), bics=bics, fullusa=fullusa,
                 state=state.in, date_run=Sys.Date(), max_date=md, lookback_days=lookback_days)
 class(ans) = "bic_seq"
 ans
}

#' report for optimal AR/MA search using BIC
#' @param x instance of bic_seq S3 class
#' @param \dots not used
#' @export
print.bic_seq = function(x, ...) {
 cat(paste("bic_seq for", x$state, "; last date used =", lubridate::as_date(x$max_date), "\n"))
 cat(" ", x$lookback_days, "day lookback.\n")
 cat(" best BIC =", min(x$bics), "for AR order", x$opt["ARord"], "MA order", x$opt["MAord"], "\n")
}

#' organize search over contiguous states for optimal statewise AR/MA using BIC 
#' @param src tibble with cumulative incidence like output of nytimes_state_data()
#' @param \dots passed to `min_bic`
#' @return instance of S3 class `min_bic_all_states`
#' @export
min_bic_all_states = function(src, ...) {
 ans = lapply(contig_states_dc(), function(x) min_bic(src, state.in=x, ...))
 names(ans) = contig_states_dc()
 class(ans) = "min_bic_all_states"
 ans
}

#' @export
print.min_bic_all_states = function(x, ...) {
 cat("min BICs for", x[[1]]$lookback_days, "day lookback from", 
    paste(lubridate::as_date(x[[1]]$max_date), "\n"))
 ars = sapply(x, function(x) x$opt["ARord"])
 mas = sapply(x, function(x) x$opt["MAord"])
 cat(" distribution of optimal autoregressive order:\n")
 print(table(ars))
 cat(" distribution of optimal moving average order:\n")
 print(table(mas))
}

#' provide vector of contiguous US state names and DC
#' @importFrom graphics abline
#' @importFrom stats coef fitted.values
#' @import dplyr
#' @import magrittr
#' @export
contig_states_dc = function() c("Alabama", "Arizona", "Arkansas", "California", "Colorado", 
"Connecticut", "Delaware", "District of Columbia", "Florida", 
"Georgia", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", 
"Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", 
"Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", 
"Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", 
"New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", 
"Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", 
"Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", 
"West Virginia", "Wisconsin", "Wyoming")

#' provide vector of contiguous US state names and DC, abbr
#' @export
contig_states_twolet = function() {
c("AL", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL", "GA", 
"ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", 
"MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", 
"ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", 
"VT", "VA", "WA", "WV", "WI", "WY")
}

# from https://worldpopulationreview.com/states/
pops = c(California = 39937489, Texas = 29472295, Florida = 21992985, 
`New York` = 19440469, Pennsylvania = 12820878, Illinois = 12659682, 
Ohio = 11747694, Georgia = 10736059, `North Carolina` = 10611862, 
Michigan = 10045029, `New Jersey` = 8936574, Virginia = 8626207, 
Washington = 7797095, Arizona = 7378494, Massachusetts = 6976597, 
Tennessee = 6897576, Indiana = 6745354, Missouri = 6169270, Maryland = 6083116, 
Wisconsin = 5851754, Colorado = 5845526, Minnesota = 5700671, 
`South Carolina` = 5210095, Alabama = 4908621, Louisiana = 4645184, 
Kentucky = 4499692, Oregon = 4301089, Oklahoma = 3954821, Connecticut = 3563077, 
Utah = 3282115, Iowa = 3179849, Nevada = 3139658, Arkansas = 3038999, 
`Puerto Rico` = 3032165, Mississippi = 2989260, Kansas = 2910357, 
`New Mexico` = 2096640, Nebraska = 1952570, Idaho = 1826156, 
`West Virginia` = 1778070, Hawaii = 1412687, `New Hampshire` = 1371246, 
Maine = 1345790, Montana = 1086759, `Rhode Island` = 1056161, 
Delaware = 982895, `South Dakota` = 903027, `North Dakota` = 761723, 
Alaska = 734002, `District of Columbia` = 720687, Vermont = 628061, 
Wyoming = 567025)


make_cumul_events = function(count, dates, 
    alpha3="USA", source="NYT", regtag=NA) {
    ans = list(count = count, dates = dates)
    attr(ans, "ProvinceState") = regtag
    attr(ans, "source") = source
    attr(ans, "alpha3") = alpha3
    attr(ans, "dtype") = "cumulative"
    class(ans) = c("cumulative_events", "covid_events")
    ans 
}

form_inc_state = function(src, regtag, max_date=NULL) {
 fullsumm = src %>% 
  dplyr::select(state,date,count) %>% group_by(date) %>% 
   summarise(count=sum(count))  # counts by date collapsed over states
 if (!is.null(max_date)) fullsumm = filter(fullsumm, date <= lubridate::as_date(max_date))
 thecum = make_cumul_events(count=fullsumm$count, dates=fullsumm$date, regtag=regtag)
 form_incident_events(thecum)
}

form_inc_county = function(src, regtag, max_date=NULL) {
 fullsumm = src %>% 
  dplyr::select(state,date,county,count) 
 if (!is.null(max_date)) fullsumm = filter(fullsumm, date <= lubridate::as_date(max_date))
 thecum = make_cumul_events(count=fullsumm$count, dates=fullsumm$date, regtag=regtag)
 form_incident_events(thecum)
}

form_inc_nation = function(src, regtag, max_date=NULL) {
 fullsumm = src %>% 
  dplyr::select(date,count) %>% group_by(date) %>% 
   summarise(count=sum(count))  # counts by date 
 if (!is.null(max_date)) fullsumm = filter(fullsumm, date <= lubridate::as_date(max_date))
 thecum = make_cumul_events(count=fullsumm$count, dates=fullsumm$date, regtag=regtag)
 form_incident_events(thecum)
}

#' Use Rob Hyndman's forecast package to estimate drift in ARIMA models
#' @import forecast
#' @importFrom lubridate as_date
#' @param src a tibble as returned by nytimes_state_data() or jhu_us_data()
#' @param state.in character(1) state name
#' @param MAorder numeric(1) order of moving average component
#' @param Difforder numeric(1) order of differencing d in ARIMA(p,d,q)
#' @param basedate character(1) used by lubridate::as_date to filter away all earlier records
#' @param lookback_days numeric(1) only uses this many days from most recent in src
#' @param ARorder order of autoregressive component
#' @param max_date a date from which to start lookback ... defaults to NULL in which
#' case the latest available date is used
#' @return instance of S3 class Arima_sars2pack
#' @examples
#' nyd = nytimes_state_data()
#' mb = min_bic(nyd, state.in="New York")
#' lkny = Arima_by_state(nyd, ARorder=mb$opt["ARord"], MAorder=mb$opt["MAord"])
#' lkny
#' plot(lkny)
#' usd = jhu_us_data()
#' lkny2 = Arima_by_state(usd, ARorder=mb$opt["ARord"], MAorder=mb$opt["MAord"])
#' lkny2
#' plot(lkny2)
#' lknyNULL = Arima_by_state(nyd, ARorder=NULL, MAorder=NULL)
#' lknyNULL
#' @export
Arima_by_state = function(src, state.in="New York", MAorder=NULL,
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=NULL, max_date=NULL) {
   cbyd = dplyr::filter(src, date >= basedate & subset==event_tag() & state==state.in) 
   ibyd = form_inc_state(cbyd, regtag=state.in, max_date=max_date)
   tc = match.call()
   if (is.null(MAorder) | is.null(ARorder) ) {
     mb = try(min_bic(src=src, state.in=state.in, basedate=basedate, lookback_days=lookback_days,
        max_date=max_date), silent=TRUE)
     MAorder = mb$opt["MAord"]
     ARorder = mb$opt["ARord"]
     }
   .Arima_inc(ibyd, state.in=state.in, MAorder=MAorder,
      Difforder=Difforder, basedate=basedate, lookback_days=lookback_days, ARorder=ARorder,
      max_date=max_date, topcall=tc)
   }

#' Use Rob Hyndman's forecast package to estimate drift in ARIMA models
#' @param ejhu a tibble as returned by `enriched_jhu_data()`
#' @param alp3 character(1) alpha3Code value for country
#' @param MAorder numeric(1) order of moving average component
#' @param Difforder numeric(1) differencing order d of ARIMA(p,d,q)
#' @param basedate character(1) used by lubridate::as_date to filter away all earlier records
#' @param lookback_days numeric(1) only uses this many days from most recent in ejhu
#' @param ARorder order of autoregressive component
#' @param max_date a date from which to start lookback ... defaults to NULL in which
#' case the latest available date is used
#' @return instance of S3 class Arima_sars2pack
#' @examples
#' ej = enriched_jhu_data()
#' lkus = Arima_nation(ej)
#' lkus
#' plot(lkus)
#' @export
Arima_nation = function(ejhu, alp3="USA", MAorder=NULL,
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=NULL, max_date=NULL ) {
   cbyd = dplyr::filter(ejhu, date >= basedate & subset==event_tag() & alpha3Code==alp3)
   ibyd = form_inc_nation(cbyd, regtag=alp3, max_date=max_date)
   if (is.null(ARorder) | is.null(MAorder)) {
     curbic = min_bic(ejhu, fullusa=TRUE)
     ARorder = curbic$opt["ARord"]
     MAorder = curbic$opt["MAord"]
     }
   tc = match.call()
   .Arima_inc(ibyd, state.in=alp3, MAorder=MAorder,
      Difforder=Difforder, basedate=basedate, lookback_days=lookback_days, ARorder=ARorder,
        max_date=max_date, topcall=tc)
   }

.Arima_inc = function(ibyd, state.in="New York", MAorder=2, 
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=0, max_date=NULL, topcall=NULL) {
   if (is.null(MAorder) | is.null(ARorder)) stop("MA/AR order inputs cannot be NULL")
   iuse = trim_from(ibyd, basedate)
   full29 = tail(ibyd$count,lookback_days)
   dates29 = tail(ibyd$date,lookback_days)
   nlb=lookback_days-1
   time = (0:nlb)/nlb

   origin = max(ibyd$date)-lookback_days+1
   time.from.origin = as.numeric(dates29-origin)
   tsfull = ts(full29, freq=1)
   Arima.full = try(Arima(tsfull, order=c(ARorder,Difforder,MAorder), include.drift=TRUE), silent=TRUE)
   if (inherits(Arima.full, "try-error")) {
     print(c(ARorder,Difforder,MAorder))
     print(Arima.full)
     ans = NA_real_
     names(ans) = "likely_AR_error"
     class(ans) = c("AR_error", "try-error")
     return(ans)
     }
   pr = fitted.values(forecast(Arima.full))
   ans = list(fit=Arima.full, pred=pr, tsfull=tsfull, dates29=dates29, time.from.origin=time.from.origin,
        call=topcall,
        state=state.in, origin=as_date(origin), MAorder=MAorder, Difforder=Difforder, ARorder=ARorder,
            max_date=max_date, Arima.inc.call=match.call())
   class(ans) = "Arima_sars2pack"
   ans
}

#' @export
print.Arima_sars2pack = function(x, ...) {
 cat("Arima_sars2pack instance for", x$state, "\n  computed", date(),  "\n")
 cat(paste("  last date used was: ", lubridate::as_date(max(x$dates29)), "\n"))
 cat("  call was: ")
 print(x$call)
 cat("  Model estimates: \n")
 summary(x$fit)
}

#' @export
plot.Arima_sars2pack = function(x, y, ...) {
# revised 8 jul to deal better with short series
    y_ = x$tsfull
    x_ = x$origin + x$time.from.origin
    yshift = floor(length(x_)/2)
    plot(x_, y_, pch = 19, xlab = "date", ylab = "incidence", 
        ...)
    lines(x$origin + x$time.from.origin, x$pred)
    if (x$Difforder == 1) {
        slo = coef(x$fit)["drift"]
        y1 = median(y_)
        y0 = y1 + slo * (-as.numeric(x$origin) - yshift)
        abline(y0, slo, lty = 2, lwd = 2)
        se = sqrt(x$fit$var.coef["drift", "drift"])
        y1b = slo * (yshift+1 + as.numeric(x$origin)) + y0
        y0p = y1b + (slo + 1.96 * se) * (-as.numeric(x$origin) - 
            yshift-1)
        y0m = y1b + (slo - 1.96 * se) * (-as.numeric(x$origin) - 
            yshift-1)
        abline(y0p, slo + 1.96 * se, lty = 3, col = "gray", lwd=2)
        abline(y0m, slo - 1.96 * se, lty = 3, col = "gray", lwd=2)
        legend("topleft", lty=c(1,2,3), lwd=c(1,2,2), 
                  legend=c("ARIMA fit", "drift component", "drift +/- 1.96sedrift"), bg="transparent")
    }
}

#' fit ARIMA model to US data dropping one state
#' @param src_us tibble for national level data like that of enriched_jhu_data()
#' @param src_st tibble for state level data like that of nytimes_state_data()
#' @param state.in character(1) state name
#' @param MAorder numeric(1) order of moving average component
#' @param Difforder numeric(1) differencing order d of ARIMA(p,d,q)
#' @param basedate character(1) used by lubridate::as_date to filter away all earlier records
#' @param lookback_days numeric(1) only uses this many days from most recent in src
#' @param ARorder order of autoregressive component
#' @param max_date a date from which to start lookback ... defaults to NULL in which
#' case the latest available date is used
#' @return instance of S3 class Arima_sars2pack
#' @note Apparent discrepancies in counts in vicinity of April 15 between
#' NYT and JHU for full US incidence lead us to employ the two sources
#' for the example.  If NYT alone is used, summing to obtain USA
#' overall, the resulting USA incidence trend is
#' extravagantly variable.
#' @examples
#' ej = enriched_jhu_data()
#' usa_full = Arima_nation(ej)
#' plot(usa_full)
#' nyd = nytimes_state_data()
#' drny = Arima_drop_state(ej, nyd)
#' drny
#' plot(drny)
#' opar = par(no.readonly=TRUE)
#' par(mar=c(4,3,2,2), mfrow=c(1,2))
#' plot(usa_full, main="all states", ylim=c(17000,38000))
#' plot(drny, main="excluding NY", ylim=c(17000,38000))
#' par(opar)
#' @export
Arima_drop_state = function(src_us, src_st, state.in="New York", MAorder=2, 
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=0, max_date=NULL) {
   curbic = min_bic(src_us, fullusa=TRUE)
   ARorder_nat = curbic$opt["ARord"]
   MAorder_nat = curbic$opt["MAord"]
   nat = Arima_nation(src_us, MAorder=MAorder_nat, Difforder=Difforder, basedate=basedate,
         lookback_days=lookback_days, ARorder=ARorder_nat, max_date=max_date)
   stbic = min_bic(src_st, fullusa=FALSE)
   ARorder_st = stbic$opt["ARord"]
   MAorder_st = stbic$opt["MAord"]
   st = Arima_by_state(src_st, state.in=state.in, MAorder=MAorder_st, Difforder=Difforder, basedate=basedate,
         lookback_days=lookback_days, ARorder=ARorder_st, max_date=max_date)
   cbyd_shim = dplyr::filter(src_st,  # shim
            date >= basedate & subset==event_tag() & state==state.in)
   ibyd_shim = form_inc_state(cbyd_shim, regtag=state.in, max_date=max_date)
   ibyd_shim$count = as.numeric(nat$tsfull)-as.numeric(st$tsfull)
   .Arima_inc(ibyd_shim, state.in=paste("excl", state.in), MAorder=MAorder_nat,
      Difforder=Difforder, basedate=basedate, lookback_days=lookback_days, ARorder=ARorder_nat,
           max_date=max_date)
   }

#' multistate exclusion prior to ARIMA
#' @param src_us tibble for national level data like that of enriched_jhu_data()
#' @param src_st tibble for state level data like that of nytimes_state_data()
#' @param states.in character() vector of state names
#' @param MAorder numeric(1) order of moving average component
#' @param Difforder numeric(1) differencing order d of ARIMA(p,d,q)
#' @param basedate character(1) used by lubridate::as_date to filter away all earlier records
#' @param lookback_days numeric(1) only uses this many days from most recent in src
#' @param ARorder order of autoregressive component
#' @param max_date character(1) or date
#' @param ARorder.nat order of autoregressive component for entire nation
#' @export
Arima_drop_states = function(src_us, src_st, states.in= c("New York", "New Jersey"), MAorder=3, 
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=0, max_date=NULL, ARorder.nat=3) {
   curbic = min_bic(src_us, fullusa=TRUE)
   ARorder_nat = curbic$opt["ARord"]
   MAorder_nat = curbic$opt["MAord"]
   nat = Arima_nation(src_us, MAorder=MAorder_nat, Difforder=Difforder, basedate=basedate,
         lookback_days=lookback_days, ARorder=ARorder_nat, max_date=max_date)
   sts = lapply(states.in, function(x) Arima_by_state(src_st, state.in=x, 
         MAorder=MAorder, Difforder=Difforder, basedate=basedate,
         lookback_days=lookback_days, ARorder=ARorder, max_date=max_date))
   names(sts) = states.in
   cbyd_shims = lapply(states.in, function(x) dplyr::filter(src_st,  # shim
            date >= basedate & subset==event_tag() & state==x))
   names(cbyd_shims) = states.in
   ibyd_shims = lapply(states.in, function(x) form_inc_state(cbyd_shims[[x]], regtag=x, max_date=max_date))
   ibyd_shim = ibyd_shims[[1]]
   ibyd_shim$count = as.numeric(nat$tsfull)-as.numeric(sts[[1]]$tsfull)
   for (i in 2:length(ibyd_shims)) 
       ibyd_shim$count = as.numeric(ibyd_shim$count)-as.numeric(sts[[i]]$tsfull)
   .Arima_inc(ibyd_shim, state.in=paste("excl", paste(states.in, collapse=", ")), MAorder=MAorder_nat,
      Difforder=Difforder, basedate=basedate, lookback_days=lookback_days, ARorder=ARorder_nat,
           max_date=max_date)
   }

#' full incidence for contiguous states
#' @inheritParams Arima_drop_state
#' @param src tibble for state level data like that of nytimes_state_data()
#' @param contig_vec character() vector of tokens for subsetting src
#' @examples
#' usd = usa_facts_data()
#' cont = Arima_contig_states(usd)
#' cont
#' plot(cont)
#' @export
Arima_contig_states = function(src, state.in="All contig", MAorder=2, 
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=0,
   contig_vec = contig_states_twolet(), max_date=NULL) {
   cbyd = dplyr::filter(src, date >= basedate & 
       subset==event_tag() & state %in% contig_vec)
   ibyd = form_inc_state(cbyd, regtag=state.in, max_date=max_date)
   if (is.null(MAorder) | is.null(ARorder)) {
     curbic = min_bic(src, fullusa=TRUE)
     ARorder = curbic$opt["ARord"]
     MAorder = curbic$opt["MAord"]
     }
   .Arima_inc(ibyd, state.in="all", MAorder=MAorder,
      Difforder=Difforder, basedate=basedate, lookback_days=lookback_days, ARorder=ARorder, max_date=max_date)
   }

#' county-level model
#' @param src a tibble as returned by nytimes_state_data() or jhu_us_data()
#' @param state.in character(1) state name
#' @param county.in character(1) state name
#' @param MAorder numeric(1) order of moving average component, defaults to NULL in which case min_bic is used to
#' find best value on grid (0:5)^2
#' @param Difforder numeric(1) order of differencing d in ARIMA(p,d,q)
#' @param basedate character(1) used by lubridate::as_date to filter away all earlier records
#' @param lookback_days numeric(1) only uses this many days from most recent in src
#' @param ARorder order of autoregressive component, defaults to NULL for optimizing choice, see MAorder
#' @param max_date a date from which to start lookback ... defaults to NULL in which
#' case the latest available date is used
#' @return instance of S3 class Arima_sars2pack
#' @examples
#' nytc = nytimes_county_data()
#' Arima_by_county(nytc)
#' @export
Arima_by_county = function(src, state.in="Massachusetts", county.in="Norfolk", MAorder=NULL,
   Difforder=1, basedate="2020-02-15", lookback_days=29, ARorder=NULL, max_date=NULL) {
   cbyd = dplyr::filter(src, date >= basedate & subset==event_tag() & state==state.in & county==county.in) 
   ibyd = form_inc_county(cbyd, regtag=paste(county.in, "/", state.in, sep="/"), max_date=max_date)
   tc = match.call()
   if (is.null(MAorder) | is.null(ARorder)) {
     mb = min_bic(src=src, state.in=state.in, county.in=county.in, basedate=basedate, lookback_days=lookback_days,
        max_date=max_date)
     MAorder = mb$opt["MAord"]
     ARorder = mb$opt["ARord"]
     }
   .Arima_inc(ibyd, state.in=state.in, MAorder=MAorder,
      Difforder=Difforder, basedate=basedate, lookback_days=lookback_days, ARorder=ARorder,
      max_date=max_date, topcall=tc)
   }
