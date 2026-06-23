## Test Cases for RTMB Gaussian Family
## Tests follow this pattern: fit with both the RTMB and legacy TMB
## backends and compare results (logLik, fixed effects, etc.)

## Example:
# remotes::install_github("glmmTMB/glmmTMB/glmmTMB", ref = "portRTMB")

#  library(glmmTMB)
# glmmTMB:::useRTMB(TRUE)
# m1 <- glmmTMB(count ~ mined + (1|site),
#                family=poisson, data=Salamanders)
# logLik(m1)
# sim <- m1$obj$simulate(complete=TRUE)

# glmmTMB:::useRTMB(FALSE)
# m1 <- glmmTMB(count ~ mined + (1|site),
#                family=poisson, data=Salamanders)
# logLik(m1)

##SLEEP STUDY DATASET

context("RTMB Gaussian backend")

skip_if_not_installed("RTMB")

data("sleepstudy", package = "lme4")

old_use_rtmb <- glmmTMB:::useRTMB()
withr::defer(glmmTMB:::useRTMB(old_use_rtmb), testthat::teardown_env())

test_that("gaussian: fixed conditional effects only", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-6)
  expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
               tolerance = 1e-6)
})

test_that("gaussian: fixed effects with offset", {
  sleepstudy$off <- 0.1 * sleepstudy$Days

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + offset(off), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + offset(off), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
})

test_that("gaussian: weighted observations (runif weights)", {
  set.seed(102)
  sleepstudy$w <- runif(nrow(sleepstudy), 0.5, 2)

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     weights = w, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
})

test_that("gaussian: weighted fixed effects (alternating weights)", {
  sleepstudy$w <- rep(c(1, 2), length.out = nrow(sleepstudy))

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     weights = w, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-6)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = 1e-6)
})

test_that("gaussian: fixed dispersion formula", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, dispformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, dispformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-6)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = 1e-6)
  expect_equal(unname(fixef(m_rtmb)$disp), unname(fixef(m_tmb)$disp),
               tolerance = 1e-6)
})

test_that("gaussian: random effects in dispersion formula", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, dispformula = ~ 1 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, dispformula = ~ 1 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = 1e-6)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = 1e-4)
})

test_that("gaussian: fixed ZI formula (hurdle, introduced zeros)", {
  set.seed(101)
  sleepstudy$Reaction[sample(nrow(sleepstudy), 5)] <- 0

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
  expect_equal(unname(fixef(m_rtmb)$zi), unname(fixef(m_tmb)$zi),
               tolerance = 1e-5)
})

test_that("gaussian: zero-inflation fixed effects match TMB backend (no induced zeros)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-6)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = 1e-6)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = 1e-5)
})

test_that("gaussian: random-only ZI formula (~0 + RE, induced zeros)", {
  set.seed(104)
  sleepstudy$Reaction[sample(nrow(sleepstudy), 5)] <- 0

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 0 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 0 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
})

test_that("gaussian: ZI intercept + random effects match TMB backend (no induced zeros)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 1 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 1 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = 1e-6)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = 1e-4)
})

test_that("gaussian: single random intercept (cond RE)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
  expect_equal(as.numeric(VarCorr(m_rtmb)$cond$Subject),
               as.numeric(VarCorr(m_tmb)$cond$Subject),
               tolerance = 1e-4)
})

test_that("gaussian: correlated random slope (us covstruct)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (Days | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (Days | Subject), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
})

test_that("gaussian: multiple random-effect terms", {
  set.seed(103)
  sleepstudy$grp2 <- factor(sample(1:5, nrow(sleepstudy), replace = TRUE))

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | grp2),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | grp2),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
})

test_that("gaussian: diag covstruct random effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + diag(Days | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + diag(Days | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = 1e-5)
})

test_that("gaussian: simulate() works under RTMB backend", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(sleepstudy))
  expect_true(all(is.finite(sim$yobs)))
})