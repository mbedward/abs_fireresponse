test_that("fesm to expert severity works for valid values", {
  fesm <- c(0, 2:5)
  expert <- c(0, 2, 4, 6, 8)
  expect_identical(fesm_to_expert_severity(fesm), expert)
})

test_that("non-integer fesm value raises error", {
  expect_error(fesm_to_expert_severity(2.5))
})

test_that("fesm 1 is treated as unburnt (0)", {
  expect_identical(fesm_to_expert_severity(1), 0)
})

test_that("out of range fesm values raise error", {
  expect_error(fesm_to_expert_severity(-1))
  expect_error(fesm_to_expert_severity(6))
})


test_that("expert to fesm severity works for valid values", {
  expert <- c(0, 2, 4, 6, 8)
  fesm <- c(0, 2, 3, 4, 5)
  expect_identical(expert_to_fesm_severity(expert), fesm)
})

test_that("non-integer expert value raises error", {
  expect_error(expert_to_fesm_severity(2.5))
})

test_that("intermediate expert values cause error by default", {
  expert <- c(1, 3, 5, 7)
  for (ex in expert) {
    expect_error(expert_to_fesm_severity(expert))
  }
})

test_that("intermediate expert values can be mapped to lower FESM values", {
  expert <- c(1, 3, 5, 7)
  fesm <- c(0, 2, 3, 4)
  expect_identical(expert_to_fesm_severity(expert, intermediate = "down"), fesm)
})

test_that("intermediate expert values can be mapped to higher FESM values", {
  expert <- c(1, 3, 5, 7)
  fesm <- c(2, 3, 4, 5)
  expect_identical(expert_to_fesm_severity(expert, intermediate = "up"), fesm)
})

test_that("out of range expert values raise error", {
  expect_error(expert_to_fesm_severity(-1))
  expect_error(expert_to_fesm_severity(9))
})



