test_that("Individual model functions", {
  data <- survival::lung
  mod <- build_cox_model(data, 'time', 'status', c('age', 'sex'))
  expect_equal(mod$nevent, 165)
  mod2 <- build_reg_model(data, 'age', 'sex')
  expect_equal(mod2$df.residual, 226)
  expect_equal(calculate_Uno_c(data, mod), 0.5848518, tolerance = 0.0001)
  expect_equal(redundancy_analysis(mod2, data)$In, c('age', 'sex'))
  expect_equal(table_predictors(data, mod, 'sex')$coef, -0.5132185, tolerance = 0.0001)
})

test_that("Blanket stats /-ments", {
  data <- survival::lung
  models_to_run <- list('OS' = list('outcome' = 'time', 'modality' = 'cox', 'event_censor' = 'status'),
                        'weight_loss' = list('outcome' = 'wt.loss', 'modality' = 'linear', 'event_censor' = NA))
  predictor_sets <- list('age' = c('age'),
                         'age_ecog' = c('age', 'ph.ecog'))
  covariates = c('sex')
  bl_stats <- blanket_statsments(data, models_to_run, predictor_sets, covariates)
  bl_redun <- blanket_redundancy_analysis(bl_stats, data)
  expect_equal(bl_stats$OS$age$nevent, 165) # if blanket_statsments works, so does blanket_stats
  expect_equal(table_blanket_statsments(data, bl_stats)['weight_loss', 'age_ecog_R^2'], 0.04990419, tolerance = 0.0001)
  expect_equal(bl_redun$OS$age$rsq1[['time']], 0.02072846 , tolerance=0.0001)
  expect_equal(table_blanket_redundancies(bl_redun)$age_ecog_redundant_vars, c('', ''))
})
