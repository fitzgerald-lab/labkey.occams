library(openclinica.occams)

test_that("connection works", {
  testthat::expect_error(connect.to.labkey("foo.txt"))
} )


test_that("demographics", {
  ocs <- connect.to.labkey()

  table = list.clinical.tables(ocs, 'di')

  di = read.demographics(ocs, table, rulesFile = 'inst/editrules/di_editrules.txt')  
  
  testthat::expect_true( is_tibble(di) & ncol(di) == 10 )
})


test_that("exposures", {
  ocs <- connect.to.labkey()
  
  ex = read.exposures(ocs, list.clinical.tables(ocs, 'ex'), rulesFile = NULL)  
  
  testthat::expect_true( (is.list(ex) && length(ex) == 2 && names(ex) %in% c('ex','history') && ncol(ex$ex) == 71) )
})

test_that('prestage', {
  ocs <- connect.to.labkey()
  
  ps = read.prestage(ocs, list.clinical.tables(ocs, 'ps'), rulesFile = NULL)  
  
  testthat::expect_true( (is.list(ps) && length(ps) == 2 && names(ps) %in% c('ps','duplicates') && ncol(ps$ps) == 65) )
})


test_that('prestage', {
  ocs <- connect.to.labkey()
  
  ps = read.prestage(ocs, list.clinical.tables(ocs, 'ps'), rulesFile = NULL)  
  
  testthat::expect_true( (is.list(ps) && length(ps) == 2 && names(ps) %in% c('ps','duplicates') && ncol(ps$ps) == 65) )
})


test_that('treatment plan', {
  ocs <- connect.to.labkey()

  tp = read.treatment.plan(ocs, list.clinical.tables(ocs, 'tp'), rulesFile = NULL)  
  
  testthat::expect_true( (is.list(tp) && length(tp) == 1 && names(tp) %in% c('tp') && ncol(tp$tp) == 13) )
})


test_that('therapy', {
  ocs <- connect.to.labkey()
  
  tr = read.therapy(ocs, list.clinical.tables(ocs, 'tr'), rulesFile = NULL)  
  
  testthat::expect_true( (is_tibble(tr) && ncol(tr) == 74) )
})


test_that('surgery', {
  ocs <- connect.to.labkey()
  
  st = read.surgery(ocs, list.clinical.tables(ocs, 'st'), rulesFile = NULL)  
  
  testthat::expect_true( (is_tibble(st) && ncol(st) == 23) )
})


test_that('pathology', {
  ocs <- connect.to.labkey()
  
  rp = read.pathology(ocs, list.clinical.tables(ocs, 'rp'), rulesFile = 'inst/editrules/rp_editrules.txt', verbose=T)  
  
  testthat::expect_true( (is_tibble(rp) && ncol(rp) == 39) )
  testthat::expect_equal(length(grep(".*\\.c", names(rp))),2)
})


test_that('followup, endpoint & recurrence', {
  ocs <- connect.to.labkey()

  fe = read.followup(ocs, tables = NULL, rulesFiles = c('inst/editrules/fe_editrules.txt', 'inst/editrules/fe1_editrules.txt', 'inst/editrules/ep_editrules.txt'))
    
  testthat::expect_true( (is_tibble(fe$fe) && ncol(fe$fe) == 13) )
  testthat::expect_equal(length(grep("Patient.Died.c", names(fe$fe))),1)
})


