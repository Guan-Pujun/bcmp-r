test_that("batch reference uses sorted levels and Python-compatible encodings", {
  ref <- bcmp:::prepare_bcmp_batch_reference(c("b", "a", "b", "c", "a"))
  expect_identical(ref$batch_levels, c("a", "b", "c"))
  expect_identical(ref$batch_idx, c(1L, 0L, 1L, 2L, 0L))
  expect_identical(ref$batch_counts, c(2L, 2L, 1L))
})

test_that("default and explicit batch coverage follow the frozen contract", {
  expect_identical(bcmp:::required_batch_coverage(3L), 2L)
  expect_identical(bcmp:::required_batch_coverage(5L), 3L)
  expect_identical(bcmp:::required_batch_coverage(5L, 4L), 4L)
  expect_error(bcmp:::required_batch_coverage(3L, 1L), "at least 2")
})

test_that("mixing rejects invalid batch configurations like Python", {
  expect_error(bcmp:::required_batch_coverage(1L), "at least two batches")
  expect_error(bcmp:::required_batch_coverage(3L, 4L), "cannot exceed")
  expect_error(bcmp:::prepare_bcmp_batch_reference(character()), "at least one batch label")
  expect_error(bcmp:::prepare_bcmp_batch_reference(c("one", "one")), "at least two batches")
})

test_that("mixing rejects invalid evaluation inputs like Python", {
  domains <- c("a", "a", "b", "b")
  batches <- c("x", "x", "y", "y")
  expect_error(
    bcmp:::evaluate_bcmp_batch_mixing_labels(domains),
    "Either batch_labels or batch_reference"
  )
  expect_error(
    bcmp:::evaluate_bcmp_batch_mixing_labels(domains, c(batches, "x")),
    "length must match"
  )
  expect_error(
    bcmp:::evaluate_bcmp_batch_mixing_labels(c("a", "b"), c("x", "y")),
    "Too few cells"
  )
  expect_error(
    bcmp:::evaluate_bcmp_batch_mixing_labels(domains, batches,
      min_cells_in_domain_per_batch = 0L),
    "positive integer"
  )
  expect_error(
    bcmp:::evaluate_bcmp_batch_mixing_labels(domains, batches,
      max_residual_cell_frac = 1),
    "finite numeric fraction"
  )
})

test_that("mixing evaluation produces the expected passing domain summary", {
  domains <- c("A", "A", "A", "A", "B", "B", "B", "C", "C")
  batches <- c("x", "x", "y", "y", "x", "y", "z", "x", "z")
  mix <- bcmp:::evaluate_bcmp_batch_mixing_labels(
    domains,
    batches,
    min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0
  )

  expect_true(mix$pass)
  expect_identical(mix$mixing_status, "pass")
  expect_identical(mix$fail_domain_count, 0L)
  expect_identical(mix$n_evaluated_domains, 3L)
  expect_identical(mix$per_domain$domain, c("A", "B", "C"))
  expect_identical(mix$per_domain$effective_batch_count, c(2L, 3L, 2L))
  expect_identical(mix$per_domain$domain_pass, c("pass", "pass", "pass"))
})

test_that("mixing residual selection follows size then first-seen order", {
  mix <- bcmp:::evaluate_bcmp_batch_mixing_labels(
    domain_labels = c(rep("A", 5L), rep("B", 3L), rep("C", 2L)),
    batch_labels = c(rep("x", 5L), rep("y", 3L), "x", "y"),
    min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0.4
  )

  expect_identical(mix$n_evaluated_domains, 2L)
  expect_identical(mix$n_residual_domains, 1L)
  expect_identical(mix$n_residual_cells, 2L)
  expect_identical(mix$per_domain$domain_pass, c("fail", "fail", "residual"))
})

test_that("mixing flags evaluated insufficiently represented domains", {
  mix <- bcmp:::evaluate_bcmp_batch_mixing_labels(
    domain_labels = c(rep("A", 4L), rep("B", 4L)),
    batch_labels = c(rep("x", 4L), rep("y", 4L)),
    min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0
  )

  expect_false(mix$pass)
  expect_identical(mix$mixing_status, "fail")
  expect_identical(mix$fail_domain_count, 2L)
  expect_identical(mix$mixing_failure_reason, "domain_mixing_fail")
})

test_that("fractional underrepresentation folds use the unchanged threshold formula", {
  labels <- rep(c("A", "B"), each = 20)
  batches <- c(rep("x", 4), rep("y", 16), rep("x", 16), rep("y", 4))
  evaluate <- function(fold) bcmp:::evaluate_bcmp_batch_mixing_labels(
    labels, batches, min_cells_in_domain_per_batch = 1L,
    max_underrepresentation_fold = fold, max_residual_cell_frac = 0)
  expect_false(evaluate(2)$pass)
  expect_true(evaluate(2.5)$pass)
  expect_true(evaluate(3.33)$pass)
  expect_false(evaluate(1)$pass)
  expect_identical(evaluate(10), evaluate(10L))
})
