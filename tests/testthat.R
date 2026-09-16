# Keep checks within two cores; child R sessions inherit these limits.
Sys.setenv(
  RCPP_PARALLEL_NUM_THREADS = "2",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

library(testthat)
library(bcmp)

test_check("bcmp")
