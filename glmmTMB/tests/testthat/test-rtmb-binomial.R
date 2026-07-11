## Test cases for the RTMB binomial family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB Binomial backend")

skip_if_not_installed("RTMB")

data("cbpp", package = "lme4")

old_use_rtmb <- glmmTMB:::useRTMB()
testthat::teardown(glmmTMB:::useRTMB(old_use_rtmb))

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(301)
binom_dat <- data.frame(
  y = rbinom(120, size = 1, prob = 0.35),
  x = rnorm(120),
  g = factor(rep(seq_len(24), each = 5))
)

test_that("binomial: binary response with fixed effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
})

test_that("binomial: grouped response with weights", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    incidence / size ~ period,
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    incidence / size ~ period,
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
})

test_that("binomial: cbind response", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(incidence, size - incidence) ~ period,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(incidence, size - incidence) ~ period,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
})

test_that("binomial: conditional random intercept", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    incidence / size ~ period + (1 | herd),
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    incidence / size ~ period + (1 | herd),
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  expect_equal(
    as.numeric(VarCorr(m_rtmb)$cond$herd),
    as.numeric(VarCorr(m_tmb)$cond$herd),
    tolerance = tol_varcorr
  )
})

test_that("binomial: cloglog link", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = binomial(link = "cloglog"),
    data = binom_dat,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = binomial(link = "cloglog"),
    data = binom_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
})

test_that("binomial: zero-inflation fixed effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
})

test_that("binomial: simulate works under RTMB backend", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    incidence / size ~ period + (1 | herd),
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(cbpp))
  expect_true(all(sim$yobs >= 0))
  expect_true(all(sim$yobs <= cbpp$size))
})
