## ----doload-------------------------------------------------------------------
suppressPackageStartupMessages({
library(sars2app)
})


## ----lkd1, cache=FALSE--------------------------------------------------------
ej = enriched_jhu_data()
usa_full = Arima_nation(ej, max_date="2020-05-20")
usa_full


## ----lkd2, cache=FALSE--------------------------------------------------------
nyd = nytimes_state_data()
drny = Arima_drop_state(ej, nyd, max_date="2020-05-20")
drny


## ----lkd, fig.width=11, fig.height=5,cache=FALSE------------------------------
par(mfrow=c(1,2), mar=c(4,3,2,2))
plot(usa_full, main="USA aggregated", ylim=c(15900,36000))
plot(drny, main="USA excluding NY", ylim=c(15900,36000) )


## ----doall,cache=FALSE--------------------------------------------------------
allst = contig_states_dc()
allarima = lapply(allst, function(x) Arima_by_state(nyd, x))
names(allarima) = allst
drifts = sapply(allarima, function(x) coef(x$fit)["drift"])
searima = function(a) sqrt(a$fit$var.coef["drift", "drift"])
se.drifts = sapply(allarima, searima)


## ----lkmeta-------------------------------------------------------------------
library(rmeta)
meta.summaries(drifts, se.drifts)

## ----dorm, fig.height=10------------------------------------------------------
o = order(drifts)
metaplot(drifts[o], se.drifts[o], labels=allst[o], cex=.7, 
  xlab="Incidence velocity (CHANGE in number of confirmed cases/day)", ylab="State")
segments(rep(-350,46), seq(-49,-4), rep(-50,46), seq(-49,-4), lty=3, col="gray")


## ----lklkaaa------------------------------------------------------------------
data(min_bic_2020_05_20)
min_bic_2020_05_20
min_bic_2020_05_20[["New York"]]

