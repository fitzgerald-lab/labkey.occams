library(openclinica.occams)

test_that("connection works", {
  ocs = connect.to("~/tmp/labkey.txt")
  expect_type(ocs, 'OCCAMSLabkey')
expect_error(connect.to("foo.txt"))

} )
