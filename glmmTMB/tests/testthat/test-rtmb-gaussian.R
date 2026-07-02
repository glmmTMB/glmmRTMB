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
testthat::teardown(glmmTMB:::useRTMB(old_use_rtmb))

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

test_that("gaussian: fixed conditional effects only", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
               tolerance = tol_fixef)
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
               tolerance = tol_logLik)
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
               tolerance = tol_logLik)
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
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian: fixed dispersion formula", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, dispformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, dispformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(unname(fixef(m_rtmb)$disp), unname(fixef(m_tmb)$disp),
               tolerance = tol_fixef)
})

test_that("gaussian: random effects in dispersion formula", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, dispformula = ~ 1 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, dispformula = ~ 1 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
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
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$zi), unname(fixef(m_tmb)$zi),
               tolerance = tol_fixef)
})

test_that("gaussian: zero-inflation fixed effects match TMB backend (no induced zeros)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
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
               tolerance = tol_logLik)
})

test_that("gaussian: ZI intercept + random effects match TMB backend (no induced zeros)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 1 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 1 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian: dZI simulation generates structural zeros", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1,
    family = gaussian,
    data = sleepstudy,
    start = list(betazi = qlogis(0.5)),
    map = list(betazi = factor(NA)),
    se = FALSE
  )

  set.seed(1001)
  sim <- m_rtmb$obj$simulate(complete = TRUE)$yobs

  expect_length(sim, nrow(sleepstudy))
  expect_true(all(is.finite(sim)))
  expect_true(any(sim == 0))
  expect_true(any(sim != 0))
})

test_that("gaussian: single random intercept (cond RE)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(as.numeric(VarCorr(m_rtmb)$cond$Subject),
               as.numeric(VarCorr(m_tmb)$cond$Subject),
               tolerance = tol_varcorr)
})

test_that("gaussian: correlated random slope (us covstruct)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (Days | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (Days | Subject), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
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
               tolerance = tol_logLik)
})

test_that("gaussian: diag covstruct random effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + diag(Days | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + diag(Days | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
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






## Salamander Dataset

data("Salamanders", package = "glmmTMB")
Salamanders$lcount <- log1p(Salamanders$count)  ## continuous response for gaussian

test_that("gaussian (Salamanders): fixed effects with binary factor predictor", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                     se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                    se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
               tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): fixed effects with multi-level factor predictor", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + spp, family = gaussian,
                     data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + spp, family = gaussian,
                    data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
               tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): single random intercept by site", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (1 | site), family = gaussian,
                     data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + (1 | site), family = gaussian,
                    data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(as.numeric(VarCorr(m_rtmb)$cond$site),
               as.numeric(VarCorr(m_tmb)$cond$site),
               tolerance = tol_varcorr)
})

test_that("gaussian (Salamanders): crossed random intercepts (site + spp)", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (1 | site) + (1 | spp),
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + (1 | site) + (1 | spp),
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian (Salamanders): nested random slope by site", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (mined | site), family = gaussian,
                     data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + (mined | site), family = gaussian,
                    data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian (Salamanders): dispersion varying by mined status", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, dispformula = ~ mined,
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, dispformula = ~ mined,
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$disp), unname(fixef(m_tmb)$disp),
               tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): zero-inflation with factor predictor (hurdle)", {
  ## lcount already has structural zeros from count == 0 observations,
  ## so no artificial zero injection needed here
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, ziformula = ~ mined,
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, ziformula = ~ mined,
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$zi), unname(fixef(m_tmb)$zi),
               tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): zero-inflation with random intercept by site", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, ziformula = ~ 1 + (1 | site),
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, ziformula = ~ 1 + (1 | site),
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): weighted observations", {
  set.seed(201)
  Salamanders$w <- runif(nrow(Salamanders), 0.5, 2)

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                     weights = w, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                    weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian (Salamanders): simulate() works under RTMB backend", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (1 | site), family = gaussian,
                     data = Salamanders, se = FALSE)
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(Salamanders))
  expect_true(all(is.finite(sim$yobs)))
})







## ChickWeight Dataset

data("ChickWeight", package = "datasets")
chick_dat <- transform(ChickWeight,
                       Chick = factor(Chick),
                       Diet = factor(Diet))
chick_dat$off <- 0.05 * chick_dat$Time
chick_dat$w <- rep(c(1, 1.5), length.out = nrow(chick_dat))
chick_dat$grp2 <- factor(as.integer(chick_dat$Chick) %% 6)

test_that("gaussian ChickWeight: fixed conditional effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                    data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                   data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: fixed effects with offset", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + offset(off), family = gaussian,
                    data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + offset(off), family = gaussian,
                   data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian ChickWeight: weighted observations", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                    data = chick_dat, weights = w, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                   data = chick_dat, weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: fixed dispersion formula", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, dispformula = ~ Time,
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, dispformula = ~ Time,
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: conditional random intercept", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian ChickWeight: conditional random slope", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + (Time | Chick),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + (Time | Chick),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian ChickWeight: diag covariance random effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + diag(Time | Chick),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + diag(Time | Chick),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian ChickWeight: multiple random-effect terms", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick) + (1 | grp2),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick) + (1 | grp2),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian ChickWeight: fixed ZI formula with induced zeros", {
  chick_zi <- chick_dat
  set.seed(201)
  chick_zi$weight[sample(nrow(chick_zi), 10)] <- 0

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ Time,
                    family = gaussian, data = chick_zi, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ Time,
                   family = gaussian, data = chick_zi, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: ZI random effects with induced zeros", {
  chick_zi <- chick_dat
  set.seed(202)
  chick_zi$weight[sample(nrow(chick_zi), 10)] <- 0

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ 1 + (1 | Chick),
                    family = gaussian, data = chick_zi, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ 1 + (1 | Chick),
                   family = gaussian, data = chick_zi, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})





## Existing Gaussian test coverage ported to RTMB/TMB comparisons

test_that("gaussian existing basics: intercept-only random effect", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ 1 + (1 | Subject), data = sleepstudy,
                    se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ 1 + (1 | Subject), data = sleepstudy,
                   se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing basics: split random intercept and slope", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                    data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                   data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing basics: double-bar random effects", {
  data("sleepstudy", package = "lme4")

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ 1 + (Days || Subject), data = sleepstudy,
                    se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ 1 + (Days || Subject), data = sleepstudy,
                   se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing basics: bar/double-bar bug model", {
  set.seed(1)
  n <- 100
  xdata <- data.frame(
    rfac1 = as.factor(sample(letters[1:10], n, replace = TRUE)),
    rfac2 = as.factor(sample(letters[1:10], n, replace = TRUE)),
    cov = rnorm(n),
    rv = rpois(n, lambda = 2)
  )

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(rv ~ cov + (1 + cov || rfac1) + (1 | rfac2),
                    family = gaussian, data = xdata, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(rv ~ cov + (1 + cov || rfac1) + (1 | rfac2),
                   family = gaussian, data = xdata, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing Anova case: fixed dispersion indicator", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject),
                    dispformula = ~ I(Days > 5), data = sleepstudy,
                    REML = FALSE, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject),
                   dispformula = ~ I(Days > 5), data = sleepstudy,
                   REML = FALSE, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing methods: Salamanders dispersion by cover", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(count ~ cover, family = gaussian,
                    dispformula = ~ cover, data = Salamanders, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(count ~ cover, family = gaussian,
                   dispformula = ~ cover, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("gaussian existing methods: mtcars random dispersion effect", {
  mtcars$cyl <- factor(mtcars$cyl)

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(mpg ~ wt, dispformula = ~ 1 + (1 | cyl),
                    data = mtcars, family = gaussian, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(mpg ~ wt, dispformula = ~ 1 + (1 | cyl),
                   data = mtcars, family = gaussian, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing varstruc: homdiag random effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + homdiag(Days | Subject),
                    data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + homdiag(Days | Subject),
                   data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing predict: polynomial fixed effects", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ poly(Days, 3), data = sleepstudy,
                    se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ poly(Days, 3), data = sleepstudy,
                   se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian existing predict: polynomial with random intercept", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ (1 | Subject) + poly(Days, 3),
                    data = sleepstudy, se = FALSE)

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ (1 | Subject) + poly(Days, 3),
                   data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})


## Testing covariance structure options

test_that("gaussian: heterogeneous compound-symmetry covariance", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + cs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + cs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("gaussian: homogeneous compound-symmetry covariance", {
  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + homcs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + homcs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("gaussian: heterogeneous Toeplitz covariance", {
  toep_data <- sleepstudy
  toep_data$Reaction <- ifelse(toep_data$Reaction > 250, 1, 0)
  toep_data$Days <- cut(
    toep_data$Days,
    breaks = c(0, 3, 6, 10),
    right = FALSE
  )

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ toep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ toep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("gaussian: homogeneous Toeplitz covariance", {
  toep_data <- sleepstudy
  toep_data$Reaction <- ifelse(toep_data$Reaction > 250, 1, 0)
  toep_data$Days <- cut(
    toep_data$Days,
    breaks = c(0, 3, 6, 10),
    right = FALSE
  )

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ homtoep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ homtoep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("gaussian: proportional covariance", {
  matrix_names <- c("(Intercept)", "Days")
  proportional_matrix <- diag(2)
  dimnames(proportional_matrix) <- list(matrix_names, matrix_names)

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + propto(Days | Subject, proportional_matrix),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + propto(Days | Subject, proportional_matrix),
    family = gaussian,
    data = sleepstudy,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("gaussian: fixed equal-to covariance", {
  matrix_names <- c("(Intercept)", "Days")
  fixed_covariance <- matrix(c(900, 5, 5, 25), 2, 2)
  dimnames(fixed_covariance) <- list(matrix_names, matrix_names)

  glmmTMB:::useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + equalto(Days | Subject, fixed_covariance),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB:::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + equalto(Days | Subject, fixed_covariance),
    family = gaussian,
    data = sleepstudy,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("RTMB uses TMB's unstructured theta ordering in dimension four", {
  L <- matrix(
    c(
      1.0, 0.0, 0.0, 0.0,
      0.2, 1.0, 0.0, 0.0,
      0.3, 0.4, 1.0, 0.0,
      0.5, 0.6, 0.7, 1.0
    ),
    nrow = 4,
    byrow = TRUE
  )
  scale <- diag(c(0.5, 1.0, 1.5, 2.0))
  Sigma <- scale %*% tcrossprod(L) %*% scale
  theta <- glmmTMB:::as.theta.vcov(Sigma)

  equalto_term <- list(
    blockSize = 4,
    blockReps = 1,
    blockNumTheta = 10,
    blockCode = structure(15, names = "equalto"),
    fullCor = 1
  )
  equalto_result <- glmmTMB:::termwise_nll(
    numeric(4), theta, equalto_term
  )

  expect_equal(equalto_result$corr, cov2cor(Sigma), tolerance = 1e-12)
  expect_equal(equalto_result$sd, sqrt(diag(Sigma)), tolerance = 1e-12)

  loglambda <- log(2.5)
  propto_term <- equalto_term
  propto_term$blockNumTheta <- 11
  propto_term$blockCode <- structure(11, names = "propto")
  propto_result <- glmmTMB:::termwise_nll(
    numeric(4), c(theta, loglambda), propto_term
  )
  propto_cov <- propto_result$corr *
    tcrossprod(propto_result$sd)

  expect_equal(propto_cov, exp(loglambda) * Sigma, tolerance = 1e-12)
})

test_that("RTMB covariance reports respect full_cor = FALSE", {
  make_term <- function(name, block_size, block_num_theta) {
    list(
      blockSize = block_size,
      blockReps = 1,
      blockNumTheta = block_num_theta,
      blockCode = structure(
        unname(glmmTMB:::.valid_covstruct[[name]]),
        names = name
      ),
      fullCor = 0
    )
  }

  cases <- list(
    us = list(n = 2, theta = c(0, 0, 0.2)),
    cs = list(n = 2, theta = c(0, 0, 0.2)),
    homcs = list(n = 2, theta = c(0, 0.2)),
    toep = list(n = 3, theta = c(0, 0, 0, 0.2, 0.1)),
    homtoep = list(n = 3, theta = c(0, 0.2, 0.1)),
    propto = list(n = 2, theta = c(0, 0, 0.2, 0))
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    term <- make_term(name, case$n, length(case$theta))
    result <- glmmTMB:::termwise_nll(
      numeric(case$n), case$theta, term
    )

    expect_identical(
      result$corr,
      matrix(NaN, 1, 1),
      info = name
    )
  }

  ## The C++ equalto implementation always reports its fixed correlation
  ## matrix, even when full_cor is false.
  equalto_term <- make_term("equalto", 2, 3)
  equalto_result <- glmmTMB:::termwise_nll(
    numeric(2), c(0, 0, 0.2), equalto_term
  )
  expect_equal(dim(equalto_result$corr), c(2, 2))
})

test_that("RTMB full_cor = FALSE reporting works through MakeADFun", {
  glmmTMB:::useRTMB(TRUE)
  fit <- glmmTMB(
    Reaction ~ Days + (Days | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE,
    control = glmmTMBControl(full_cor = FALSE)
  )

  expect_identical(fit$obj$report()$corr[[1]], matrix(NaN, 1, 1))
})
