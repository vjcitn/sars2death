library(sars2death)
library(testthat)
library(lubridate)

Arima_state_result_names = c("fit", "pred", "tsfull", "dates29", "time.from.origin", "call", 
"state", "origin", "MAorder", "Difforder", "ARorder", "max_date")

context("basic data and accumulation tests")



test_that("nytimes state data columns are as expected", {
nyd <<- nytimes_state_data()
expect_true(all(names(nyd)==c("date", "state", "fips", "count", "subset")))
})

test_that("cumulative events for state work", {
 cc = cumulative_events_nyt_state(nyd)
 expect_true(all(class(cc)==c("cumulative_events", "covid_events")))
 expect_true(all(names(cc) == c("count", "dates")))
})

test_that("usa_facts data are as expected", {
us = usa_facts_data()
expect_true(all(c("fips", "county", "state", "subset", "date", "count") %in% names(us)))
})

context("basic Arima runs")

n1 = Arima_by_state(nyd)

test_that("returned object has expected fields", {
 expect_true(all(Arima_state_result_names %in% names(n1)))
})

test_that("latest default date close to current date", {
 expect_true((as_date(Sys.Date()) - max(n1$dates29)) < 3)
})

context("date control for lookback works")

n2 = Arima_by_state(nyd, max_date="2020-05-15")

test_that("max_date setting works", {
 expect_true((as_date("2020-05-15") == as_date(max(n2$dates29))))
})


context("test server.R code for app")


 basedate = "2020-02-15" # data prior to this date are dropped completely
 lookback_days = 29 # from present date
 if (!exists(".nyd.global")) .nyd.global <<- nytimes_state_data() # cumulative
 if (!exists(".jhu.global")) .jhu.global <<- enriched_jhu_data() # cumulative
 allst = sort(unique(.nyd.global$state))
 data(list="min_bic_2020_11_26_death", package="sars2death")

known_Arima_sars2pack_components = c("fit", "pred", "tsfull", "dates29", "time.from.origin", "call", 
"state", "origin", "MAorder", "Difforder", "ARorder", "max_date"
)

test_that("USA fit gives expected class", {
   input = list()
   input$source = "fullusa"
   input$excl = "no"
   input$Difforder=1
   input$MAorder = 3
   input$ARorder = 0
   if (input$source == "fullusa" & input$excl == "no") 
       curfit = Arima_nation(.jhu.global, Difforder=input$Difforder, MAorder=input$MAorder, 
            ARorder=input$ARorder, max_date=input$maxdate)
   expect_true(class(curfit) == "Arima_sars2pack")
   expect_true(all(known_Arima_sars2pack_components %in% names(curfit)))
   })

test_that("drop state succeeds", {
   input = list()
   input$source = "fullusa"
   input$excl = "New York"
   input$Difforder=1
   input$MAorder = 3
   input$ARorder = 0
   if (input$source == "fullusa" & input$excl != "no") 
        curfit = Arima_drop_state(.jhu.global, .nyd.global, state.in=input$excl, 
         Difforder=input$Difforder, MAorder=input$MAorder, 
         ARorder=input$ARorder, max_date=input$maxdate)
   expect_true(class(curfit) == "Arima_sars2pack")
})

test_that("by state succeeds", {
   input = list()
   input$source = "New York"
   input$Difforder=1
   input$MAorder = 3
   input$ARorder = 0
   curfit = Arima_by_state(.nyd.global, state.in=input$source, 
         Difforder=input$Difforder, MAorder=input$MAorder, 
         ARorder=input$ARorder, max_date=input$maxdate)
   expect_true(class(curfit) == "Arima_sars2pack")
})

test_that("run_meta succeeds", {
   data(list="min_bic_2020_11_26_death", package="sars2death")
   input = list()
   input$Difforder=1
   mchk = run_meta(.nyd.global, opt_parms=min_bic_2020_11_26_death, Difforder=input$Difforder, 
             max_date=input$maxdate)  # note that AR/MA parms from opt_parms
   expect_true(all(c("drifts", "se.drifts") %in% names(mchk)))
  })

#test_that("min_bic gives expected result", {
#  def_515 = min_bic(nyd, max_date="2020-05-15")
#  expect_true(all(def_515$opt == c("ARord"=2, "MAord"=1)))
#  expect_true(all(names(def_515$opt) == names(c("ARord"=2, "MAord"=1))))
#})

test_that("est_Rt succeeds", {
  ny = Arima_by_state(nyd, state.in="New York", max_date="2020-05-15")
  ee = est_Rt(ny)
  expect_true(abs(ee$R[22,"Mean(R)"]-0.61866)<.001)
})



#   else if (input$source != "fullusa") curfit = Arima_by_state(.nyd.global, state.in=input$source, Difforder=input$Difforder, MAorder=input$MAorder, ARorder=input$ARorder, max_date=input$maxdate)
#   validate(need(!inherits(curfit, "try-error"), "please alter AR or MA order"))
#   list(fit=curfit, pred=fitted.values(forecast(curfit$fit)), tsfull=curfit$tsfull, dates29=curfit$dates29)
#   })
#  output$traj = renderPlot({
#   ans = dofit()
#   validate(need(!inherits(ans$fit, "try-error"), "please alter AR order"))
#   plot(ans$fit)
#   })
#  output$rept = renderPrint({ 
#    ans = dofit()
#   validate(need(!inherits(ans$fit, "try-error"), "please alter AR order"))
#    ans$fit
#   })
#  dometa = reactive({
#    run_meta(.nyd.global, opt_parms=min_bic_2020_05_20, Difforder=input$Difforder, 
#            max_date=input$maxdate)  # note that AR/MA parms from opt_parms
#  })
#  output$meta.rept = renderPrint({ 
#    m1 = dometa()
#    summ1 = rmeta::meta.summaries(m1$drifts, m1$se.drifts)
#    names(m1$drifts) = gsub(".drift", "", names(m1$drifts))
#    nyind = which(names(m1$drifts) %in% c("New York", "New Jersey"))
#    summ2 = rmeta::meta.summaries(m1$drifts[-nyind], m1$se.drifts[-nyind])
#    list(overall=summ1, exclNYNJ=summ2)
#   })
#  output$metaplot = renderPlot({
#    m1 = dometa()
#    names(m1$drifts) = gsub(".drift", "", names(m1$drifts))
#    o = order(m1$drifts)
#    rmeta::metaplot(m1$drifts[o], m1$se.drifts[o], labels=names(m1$drifts)[o], cex=.7, 
#      xlab="Infection velocity (CHANGE in number of confirmed cases/day)", ylab="State")
#    segments(rep(-350,46), seq(-49,-4), rep(-50,46), seq(-49,-4), lty=3, col="gray")
#   })
