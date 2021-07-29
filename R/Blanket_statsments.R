# Wrapper methods to build common models quickly with a set of static covariables
# J. Peter Marquardt, 2021-07-27


#' Build a cox model
#'
#' Build a Cox proportional hazards model from data and meta-parameters
#'
#' @param df data.frame containing the data set
#' @param event_time character denoting column with event time
#' @param event_censor character denoting column specifying events/censoring
#' @param predictors character vector denoting columns with independent variables of interest
#' @param covariates character vector denoting columns with independent variables not of interest.
#' Covariates are mathematically identical to predictors but will be ignored in reporting
#' @param verbose logical. TRUE activates printout messages
#'
#' @return A Cox proportional hazards model
#'
#' @importFrom survival Surv coxph
#' @export
#'
#' @examples
#' data <- survival::lung
#' mod <- build_cox_model(data, 'time', 'status', c('age', 'sex'))
#'
#' @author J. Peter Marquardt
build_cox_model <- function(df, event_time, event_censor, predictors, covariates=c(), verbose = FALSE){

  variables <- append(covariates, predictors) # merging all independent variables into one vector

  # building the actual formula
  cox_form <- .build_model_formula(outcome = event_time, predictors = variables, censor_event = event_censor)
  cox_model <- survival::coxph(formula = cox_form, data = df)

  if (verbose){print(table_predictors(df, cox_model, predictors))}

  return(cox_model)

}


#' Build a generic regression model model
#'
#' Build a generic regression model from data and meta-parameters.
#' Currently only available for linear and logistic types.
#'
#' @param df data.frame containing the data set
#' @param outcome character denoting column with the outcome of interest
#' @param predictors character vector denoting columns with independent variables of interest
#' @param covariates character vector denoting columns with independent variables not of interest.
#' Covariates are mathematically identical to predictors but will be ignored in reporting
#' @param modality character designating type. Currently limited to 'linear' and 'logistic'.
#' @param verbose logical. TRUE activates printout messages
#'
#' @return A regression model of linear or logistic type
#'
#' @examples
#' mod <- build_reg_model(data.frame('outcome' = c(1,2), 'pred' = c(3,4)), 'outcome', c('pred'))
#'
#' @importFrom basecamb build_model_formula
#' @importFrom stats lm glm
#' @export
#'
#' @author J. Peter Marquardt
build_reg_model <- function(df, outcome, predictors, covariates=c(), modality='linear', verbose = FALSE){

  stopifnot(modality %in% c('linear', 'logistic'))  # check that a correct modality is specified

  variables <- append(covariates, predictors) # merging all into one vector

  # building the actual formula
  reg_form <- basecamb::build_model_formula(outcome = outcome, predictors = variables)

  # building and returning the linear of logistic regression model
  if (modality == 'linear') {
    lin_model <- stats::lm(formula = reg_form, data = df)
    if (verbose){print(table_predictors(df, lin_model, predictors))}
    return(lin_model)
  }

  if (modality == 'logistic') {
    log_model <- glm(formula = reg_form, data = df, family = 'binomial')
    if (verbose){print(table_predictors(df, log_model, predictors))}
    return(log_model)
  }

}



#' Calculate Uno's C for a given model.
#'
#' Calculate Uno's concordance statistic for any model.
#' CAVE: If you want to evaluate a model trained on a different dataset, df should be limited to the test set.
#'
#' @param df data.frame containing the data set. If evaluating independently, use the test set.
#' @param model statistical model of type coxph to be evaluated.
#' @param verbose logical. TRUE activates printout messages.
#'
#' @return double AUC value for the evaluated model on the specified data set.
#'
#' @examples
#' data <- survival::lung
#' cancer_mod <- survival::coxph(survival::Surv(time, status)~age, data = data)
#' calculate_Uno_c(data, cancer_mod)
#'
#' @importFrom survAUC UnoC
#' @importFrom survival Surv
#' @importFrom stats predict
#' @export
#'
#' @author J. Peter Marquardt
calculate_Uno_c <- function(df, model, verbose = FALSE){

  # creating a dataframe without NAs (survAUC implementation can't handle this)
  pred_vars <- names(model[['assign']])  # getting the variables used in the model
  df_no_na <- df
  for (coln in pred_vars){df_no_na <- df_no_na[!is.na(df_no_na[coln]),]}

  # extracting Survival model and predictor variables used in the model training
  Surv_model_old <- model$y

  # building Surv() model for the test dataframe
  form_string <- paste(deparse(model$formula), collapse = '')  # extracting formula as string
  surv_string <- strsplit(strsplit(form_string, 'Surv(', fixed=TRUE)[[1]][[2]], ')')[[1]][1] # Extracting insides of the Surv() part of it
  surv_params <- strsplit(gsub(" ", "", surv_string, fixed = TRUE), ',')  # extracting the parameters
  surv_time <- surv_params[[1]][1]  # assigning time variable name
  surv_cens <- surv_params[[1]][2]  # assigning cens variable name

  # ensuring that surv time and censor variable contain no missing values
  for (coln in c(surv_time, surv_cens)){df_no_na <- df_no_na[!is.na(df_no_na[coln]),]}

  Surv_model_new <- survival::Surv(df_no_na[[surv_time]], df_no_na[[surv_cens]])  # building the actual Surv()

  # using the model to predict outcomes on original df
  pred_scores <- stats::predict(model, newdata = df_no_na)

  # Calculating the actual concordance statistic
  c_stat <- survAUC::UnoC(Surv.rsp = Surv_model_old,
                 Surv.rsp.new = Surv_model_new,
                 lpnew = pred_scores)

  # returning
  if(verbose){print(c_stat)}
  return(c_stat)
}



#' Generic wrapper method to calculate C-statistics
#'
#' Calculate concordance statistics for a list of statistical models on the same data set
#'
#' @param df data.frame containing the data set. If evaluating independently, use the test set.
#' @param model_list list of statistical models of type lm, glm or coxph to be evaluated.
#' @param modality character specifying model type. Currently accepts 'linear', 'logistic', and 'cox'
#' @param verbose logical. TRUE activates printout messages.
#'
#' @return list of doubles with the AUC values for the evaluated models on the specified data set.
#'
#' @importFrom DescTools Cstat
#' @importFrom stats predict model.response model.frame
#' @export
#'
#' @author J. Peter Marquardt
blanket_c_statistic <- function(df, model_list, modality = 'logistic', verbose = FALSE) {

  stopifnot(modality %in% c('linear', 'logistic', 'cox'))  # check that a correct modality is specified
  cstat_list <- vector(mode = "list", length = 0)  # initialising list of objects we intend to output

  for (model_name in names(model_list)) {  # iterate over models in list

    if(modality == 'cox'){
      # calculate c-statistic according to Uno
      cstat_list[[model_name]] <- calculate_Uno_c(df, model_list[[model_name]], verbose = FALSE)
    }

    else{
      # calculate standard c-statistic
      cstat_list[[model_name]] <- DescTools::Cstat(x = stats::predict(model_list[[model_name]], method="response"),
                                        resp = stats::model.response(stats::model.frame(model_list[[model_name]])))
    }
  }

  # output
  if(verbose){print(cstat_list)}
  return(cstat_list)
}



#' Table model predictor performance
#'
#' Extract coefficients and p-values only for regression models and table them
#'
#' @param df data.frame containing the data set. If evaluating independently, use the test set.
#' @param model statistical model to be evaluated.
#' @param predictors vector of characters designating columns of interest. Non-specified independent variables will not be included.
#'
#' @return data.frame with coefficients and p-values for predictor variables
#'
#' @examples
#' data <- survival::lung
#' mod <- build_reg_model(data, 'age', 'sex')
#' tbl <- table_predictors(data, mod, 'sex')
#'
#' @export
#'
#' @author J. Peter Marquardt
table_predictors <- function(df, model, predictors) {

  # extract coefficients as df
  coeff_df <- as.data.frame(summary(model)$coefficient)

  # calculate how many rows we need to include
  num_predictor_rows <- 0
  for (pred in predictors) {

    if (is.factor(df[[pred]])) {  # double index brackets are essential here!
      # we'll have one row per factor level, minus one reference level
      num_predictor_rows <- num_predictor_rows + nlevels(df[[pred]]) -1 # nrow(unique(df[pred])) - 1
    }
    else {  # just one row
      num_predictor_rows <- num_predictor_rows + 1
    }
  }

  # return just the rows with variables of interest (last n)
  return(coeff_df[(nrow(coeff_df) - num_predictor_rows + 1) : nrow(coeff_df),])

}


#' Run multiple slightly different models of same type
#'
#' Run the same model (type, outcome, and covariates) with different sets of predictors
#'
#' @param df data.frame containing the data set.
#' @param outcome character designating the column with the outcome of interest
#' @param predictor_sets named list or character vectors containing columns with predictors
#' @param covariates vector of characters denoting columns with covariables
#' @param modality character denoting model type. Currently limited to 'linear', 'logistic', and 'cox'
#' @param event_censor character denoting column with censor event. For coxph models only
#' @param verbose logical. TRUE activates printout messages.
#'
#' @return named list of models
#'
#' @examples
#' data <- survival::lung
#' outcome <- 'time'
#' predictor_sets <- list('age' = c('age'),'age_ecog' = c('age', 'ph.ecog'))
#' covariates = c('sex')
#' modality <- 'cox'
#' event_censor <- 'status'
#' bl_stats <- blanket_stats(data, outcome, predictor_sets, covariates, modality, event_censor)
#'
#' @export
#'
#' @author J. Peter Marquardt
blanket_stats <- function(df, outcome, predictor_sets, covariates=c(), modality = 'linear', event_censor = NA, verbose = FALSE) {

  stopifnot(modality %in% c('linear', 'logistic', 'cox'))  # check that a correct modality is specified
  model_list <- vector(mode = "list", length = 0)  # initialising list of objects we intend to output

  if (modality == 'cox') {  # we're building cox regression models

    # iterate over sets of predictors to apply to the otherwise same model
    for (model_name in names(predictor_sets)) {

      pset <- predictor_sets[[model_name]]  # ;)

      # running the model for one specific set of predictors and attaching it to the list to return
      model_list[[model_name]] <- build_cox_model(predictors = pset,
                                                  df = df,
                                                  event_time = outcome,
                                                  event_censor = event_censor,
                                                  covariates = covariates,
                                                  verbose = verbose)
    }
  }


  else {  # we're building linear/logistic regression models

    # iterate over sets of predictors to apply to the otherwise same model
    for (model_name in names(predictor_sets)) {

      pset <- predictor_sets[[model_name]]  # #MITpuns

      # running the model for one specific set of predictors and attaching it to the list to return
      model_list[[model_name]] <- build_reg_model(predictors = pset,
                                                  df = df,
                                                  outcome = outcome,
                                                  covariates = covariates,
                                                  modality = modality,
                                                  verbose = verbose)
    }
  }


  # return list of models
  return(model_list)
}


#' Run multiple different models with different sets of predictors
#'
#' Wraps blanket_stats. Run a list of models with different modalities/outcomes for a list of different predictor sets with the same covariables.
#'
#' @param df data.frame containing the data set.
#' @param models_to_run either a named list or data.frame type, with every entry/row having the keys/columns outcome, modality, and event_censor
#' @param predictor_sets named list of lists containing the set of predictors. See blanket_stats for details
#' @param covariates vector of characters denoting columns with covariables
#' @param verbose logical. TRUE activates printout messages.
#'
#' @return named list of named lists of models
#'
#' @examples
#' data <- survival::lung
#' models_to_run <- list('OS' = list(
#' 'outcome' = 'time', 'modality' = 'cox', 'event_censor' = 'status'),
#' 'weight_loss' = list('outcome' = 'wt.loss', 'modality' = 'linear', 'event_censor' = NA))
#' predictor_sets <- list('age' = c('age'),'age_ecog' = c('age', 'ph.ecog'))
#' covariates = c('sex')
#' bl_stats <- blanket_statsments(data, models_to_run, predictor_sets, covariates)
#'
#' @export
#'
#' @author J. Peter Marquardt
blanket_statsments <- function(df, models_to_run, predictor_sets, covariates=c(), verbose = FALSE) {

  # if models_to_run is passed as a data frame, we need to convert it into a list
  if('data.frame' %in% class(models_to_run)){
    # initialising
    models_to_run_df <- models_to_run
    models_to_run <- vector(mode = "list", length = 0)
    for (r in 1:nrow(models_to_run_df)){
      model <- vector(mode = "list", length = 0)
      model[['outcome']] <- models_to_run_df$outcome[r]
      model[['modality']] <- models_to_run_df$modality[r]
      model[['event_censor']] <- models_to_run_df$event_censor[r]
      models_to_run[[models_to_run_df$model[r]]] <- model
    }
    rm(models_to_run_df)
  }

  all_models <- vector(mode = "list", length = 0)  # initialising a list to store the lists of models we'll run

  # running the models
  for (model_name in names(models_to_run)){
    model <- models_to_run[[model_name]]
    if(verbose){
      print(paste0('Models: ', model_name))
      print(model)
      }
    models <- blanket_stats(df, outcome = model$outcome,
                            predictor_sets=predictor_sets,
                            covariates=covariates,
                            modality = model$modality,
                            event_censor = model$event_censor,
                            verbose = verbose)
    all_models[[model_name]] <- models
  }

  return(all_models)
}



# function to compile key results of an array of statistical models in a data frame
# blanket_statsment_models is a list of lists of statistical models as returned by blanket_statsments()
# CAVE: The listed p-value only refers to the p-value of the last parameter entered into the statistical model
#' Table results of multiple different models with different sets of predictors
#'
#' Wraps blanket_stats. Run a list of models with different modalities/outcomes for a list of different predictor sets with the same covariables.
#'
#' @param df data.frame containing the data set.
#' @param blanket_statsment_models list of models produced by blanket_statsments()
#'
#' @return data.frame with tabled results
#'
#' @examples
#' data <- survival::lung
#' models_to_run <- list('OS' = list(
#' 'outcome' = 'time', 'modality' = 'cox', 'event_censor' = 'status'),
#' 'weight_loss' = list('outcome' = 'wt.loss', 'modality' = 'linear', 'event_censor' = NA))
#' predictor_sets <- list('age' = c('age'),'age_ecog' = c('age', 'ph.ecog'))
#' covariates = c('sex')
#' bl_stats <- blanket_statsments(data, models_to_run, predictor_sets, covariates)
#' tbl <- table_blanket_statsments(data, bl_stats)
#'
#' @seealso [blanket_statsments()] for models and [table_predictors()] for tabling results
#'
#' @importFrom utils tail
#' @importFrom stats terms nobs
#' @export
#'
#' @author J. Peter Marquardt
table_blanket_statsments <- function(df, blanket_statsment_models){

  results_df <- data.frame()  # initialising the df to be returned

  # filling the dataframe row by row, with each row representing results of models with shared modality and outcome, but  different predictor sets
  for (model_name in names(blanket_statsment_models)){

    models <- blanket_statsment_models[[model_name]]  # extracting a single list of models for one row

    # seperate paths depending on model type
    if (class(models[[1]])[1] == 'coxph'){
      c_stats_list <- blanket_c_statistic(df, models, modality = 'cox')

      for (pred in names(models)){
        results_df[model_name, paste0(pred, '_n')] <- models[[pred]]$n
        results_df[model_name, paste0(pred, '_C')] <- c_stats_list[[pred]]
        results_df[model_name, paste0(pred, '_R^2')] <- summary(models[[pred]])$rsq[['rsq']]
        # this slightly convoluted call extracts the p-value for the last parameter entered into the statistical model
        results_df[model_name, paste0(pred, '_p')] <- utils::tail(table_predictors(df, models[[pred]],
                                                                       labels(stats::terms(models[[pred]]))[length(labels(terms(models[[pred]])))])[['Pr(>|z|)']], 1)
      }
    }


    else if (class(models[[1]])[1] == 'lm'){
      c_stats_list <- blanket_c_statistic(df, models, modality = 'linear')

      for (pred in names(models)){
        results_df[model_name, paste0(pred, '_n')] <- stats::nobs(models[[pred]])
        results_df[model_name, paste0(pred, '_C')] <- c_stats_list[[pred]]
        results_df[model_name, paste0(pred, '_R^2')] <- summary(models[[pred]])$r.squared
        results_df[model_name, paste0(pred, '_p')] <- tail(table_predictors(df, models[[pred]],
                                                                       labels(terms(models[[pred]]))[length(labels(terms(models[[pred]])))])[['Pr(>|t|)']], 1)
      }
    }


    # CAVE: Only works because we don't use glm for anything other than logistic regression
    else if (class(models[[1]])[1] == 'glm'){
      c_stats_list <- blanket_c_statistic(df, models, modality = 'logistic')

      for (pred in names(models)){
        results_df[model_name, paste0(pred, '_n')] <- nobs(models[[pred]])
        results_df[model_name, paste0(pred, '_C')] <- c_stats_list[[pred]]
        results_df[model_name, paste0(pred, '_R^2')] <- with(summary(models[[pred]]), 1 - deviance/null.deviance)
        results_df[model_name, paste0(pred, '_p')] <- tail(table_predictors(df, models[[pred]],
                                                                       labels(terms(models[[pred]]))[length(labels(terms(models[[pred]])))])[['Pr(>|z|)']], 1)
      }
    }
  }

  return(results_df)
}



#' Redundancy analysis
#'
#' Perform a redundancy analysis on an existing model
#'
#' @param model a statistical regression model of class linear, logistic or coxph
#' @param data data.frame used to create the model
#' @param r2_threshold float threshold value to consider a parameter redundant
#' @param nk number of knots in splicing
#'
#' @return an object of class "redun"
#'
#' @examples
#' data <- survival::lung
#' mod <- build_reg_model(data, 'age', c('sex'))
#' redundancy_analysis(mod, data)
#'
#' @importFrom Hmisc redun
#' @importFrom stats as.formula
#' @export
#'
#' @author J. Peter Marquardt
redundancy_analysis <- function(model, data, r2_threshold=0.9, nk=0){
  # seperate paths depending on model type
  if ('coxph' %in% class(model)){
    outcome <- as.character(model$formula[[2]][2])
    predictors <- paste(as.character(model$formula[[3]])[-1], collapse = ' + ')
    form <- stats::as.formula(paste(outcome, predictors, sep='~'))
    return(Hmisc::redun(form, data = data, r2=r2_threshold, nk=nk))
  }

  else if ('lm' %in% class(model) | 'glm' %in% class(model)){
    return(Hmisc::redun(model$terms, data = data, r2=r2_threshold, nk=nk))
  }

  stop('The model needs to be a regression of type coxph, lm or glm.')
}


#' Blanket redundancy analysis
#'
#' Perform a blanket redundancy analysis on a list of existing models
#'
#' @param model_list a list of statistical regression model of class linear, logistic or coxph
#' @param data data.frame used to create the models
#' @param r2_threshold float threshold value to consider a parameter redundant
#' @param nk number of knots in splicing
#' @param verbose ctivate printouts of key findings
#'
#' @return an list of objects of class "redun"
#'
#' @examples
#' data <- survival::lung
#' models_to_run <- list(
#' 'OS' = list('outcome' = 'time', 'modality' = 'cox', 'event_censor' = 'status'),
#' 'weight_loss' = list('outcome' = 'wt.loss', 'modality' = 'linear', 'event_censor' = NA))
#' predictor_sets <- list('age' = c('age'), 'age_ecog' = c('age', 'ph.ecog'))
#' covariates = c('sex')
#' bl_stats <- blanket_statsments(data, models_to_run, predictor_sets, covariates)
#' blanket_redundancy_analysis(bl_stats, data)
#'
#' @seealso [blanket_stats()]
#'
#' @export
#'
#' @author J. Peter Marquardt
blanket_redundancy_analysis <- function(model_list, data, r2_threshold=0.9, nk=0, verbose=FALSE){

  redundancy_list <- list()
  # loop over list of list of models to create an analogous list of lists of redun objects
  for (outcome in names(model_list)) {
    redundancy_list[[outcome]] <- list()
    for (predictors in names(model_list[[outcome]])){
      redundancy_list[[outcome]][[predictors]] <- redundancy_analysis(model_list[[outcome]][[predictors]], data=data, r2_threshold=r2_threshold, nk=nk)
    }
  }

  return(redundancy_list)
}


#' Table results of blanket redundancy analysis
#'
#' Table results of a blanket redundancy analysis on a list of existing models
#'
#' @param blanket_redundancies list of lists of redun objects generated by blanket_redundancy_analysis()
#' @param digits integer number of decimals to include
#'
#' @return a data.frame tabling the key results
#'
#' @examples
#' data <- survival::lung
#' models_to_run <- list(
#' 'OS' = list('outcome' = 'time', 'modality' = 'cox', 'event_censor' = 'status'),
#' 'weight_loss' = list('outcome' = 'wt.loss', 'modality' = 'linear', 'event_censor' = NA))
#' predictor_sets <- list('age' = c('age'), 'age_ecog' = c('age', 'ph.ecog'))
#' covariates = c('sex')
#' bl_stats <- blanket_statsments(data, models_to_run, predictor_sets, covariates)
#' bl_redun <- blanket_redundancy_analysis(bl_stats, data)
#' table_blanket_redundancies(bl_redun)
#'
#' @seealso [table_predictors()], [blanket_redundancy_analysis()]
#'
#' @export
#'
#' @author J. Peter Marquardt
table_blanket_redundancies <- function(blanket_redundancies, digits=2){

  results_df <- data.frame()
  # filling the dataframe row by row, with each row representing results of models with shared modality and outcome, but  different predictor sets
  for (outcome in names(blanket_redundancies)) {
    for (predictors in names(blanket_redundancies[[outcome]])){
      red_vars <- names(blanket_redundancies[[outcome]][[predictors]][['rsquared']])
      red_r2 <- round(blanket_redundancies[[outcome]][[predictors]][['rsquared']], digits)
      results_df[outcome, paste0(predictors, '_redundant_vars')] <- paste(red_vars, collapse = '; ')
      results_df[outcome, paste0(predictors, '_r2_vals')] <-paste(red_r2, collapse = '; ')
    }
  }

  return(results_df)
}


#' Build formula for statistical models
#'
#' Build formula used in statistical models from vectors of strings. Copied from basecamb package to avoid dependency
#'
#' @param outcome character denoting the column with the outcome.
#' @param predictors vector of characters denoting the columns with the
#'   predictors.
#' @param censor_event character denoting the column with the censoring event,
#'   for use in Survival-type models.
#'
#' @return formula for use in statistical models
#'
#' @source \link[basecamb]{build_model_formula}
#'
#' @importFrom assertive.types assert_is_character
#' @importFrom survival Surv
#'
#' @author J. Peter Marquardt
.build_model_formula <- function(outcome, predictors, censor_event=NULL) {

  assertive.types::assert_is_character(outcome)
  assertive.types::assert_is_character(predictors)

  if(is.null(censor_event)) {  # standard formula
    frml <- as.formula(paste(outcome,
                             ' ~ ',
                             paste(predictors, collapse = ' + '),
                             sep = ''))
  }

  else {  # Survival-type formula
    assertive.types::assert_is_character(censor_event)
    frml <- as.formula(paste('Surv(',
                             outcome,
                             ', ',
                             censor_event,
                             ')~',
                             paste(predictors, collapse = ' + '),
                             sep = ''))
  }

  return(frml)
}
