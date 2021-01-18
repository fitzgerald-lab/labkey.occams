library(openclinica.occams)

test_that("connection works", {
  expect_error(connect.to.labkey("foo.txt"))
} )
