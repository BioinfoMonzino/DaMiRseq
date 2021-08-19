#' @title Select the best classification model
#'
#' @description This function selects the bestmodels learned by the
#' \link{DaMiR.EnsL_Train} and a tested by
#' \link{DaMiR.EnsL_Test}.
#'
#' @param df A data frame of performance metrics. At least two columns
#'  representing a specific classification metrics (e.g., Accuracy) and
#'  the number of predictors must be provided. Additionally,
#'  other classification metrics (e.g., MCC, Sensitivity,
#'  Specificity,PPV, NPV, AUC, ...) can be appended (from the third
#'  column onwards) and used for the evaluation, by correctly setting
#'  the 'metric.idx' parameter.
#' @param type.sel The method to select the best models. Only
#'  "mode","median" and "greater" values are allowed.
#'  For a specific classification metrics, "mode" selects all models
#'  whose score is the mode of all scores; "median" selects all
#'  models whose score is the median of all scores; and, "greater"
#'  selects all models whose score is greater than the value
#'  specified in "th.sel". Default: "mode".
#' @param th.sel Threshold for the evaluation of the performance
#'  when "type.sel" is equal to "greater". Default: 0.85
#' @param npred.sel The method to select the best model. Only
#'  "min" and "rnd" values are allowed. Taking into account the subset
#'  of models found by 'type.sel', this parameter selects one single
#'  model with the minimum number of predictors ("min") or randomly
#'  ("rnd"). Default: "min".
#'
#' @param metric.idx The index of the 'df' column (i.e., classification
#' metrics) to be considered for the models evaluation. Default: 1.
#' @param npred.idx The index of the 'df' column representing the
#' number of predictors. Default: 2.
#'
#' @return The index of df (row), representing the model selected and a
#' bubble chart
#'
#' @details
#' This function finds the best model, taking into account specific
#' classification metrics.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' set.seed(1)
#' @export
#'
#'
DaMiR.ModelSelect <- function(df,
                              type.sel=c("mode",
                                         "median",
                                         "greater"),
                              th.sel=0.85,
                              npred.sel=c("min",
                                          "rnd"),
                              metric.idx=1,
                              npred.idx=2){

  # check missing arguments
  if (missing(df))
    stop("'data' argument must be provided")
  if (missing(type.sel)){
    type.sel <- "mode"
  }
  if (missing(npred.sel)){
    npred.sel <- "min"
  }
  # check the type of argument
  if(!(is.data.frame(df)))
    stop("'df' must be a data frame")
  if(!(is.numeric(metric.idx)))
    stop("'metric.idx' must be numeric")
  if(!(is.numeric(th.sel)))
    stop("'th.sel' must be numeric")

  # specific checks

  if (ncol(df) < 2)
    stop("df must have at least 2 columns:
         the first must represent a classification metrics and
         the second th number of predictors")
  if (th.sel > 1 | th.sel < 0)
    stop("'th.sel' must be between 0 and 1")
  if (metric.idx > ncol(df) | metric.idx <= 0)
    stop("'metric.idx' must be at least equal to ncol(df)")
  if (metric.idx <= 0)
    stop("'metric.idx' must be positive")
  if (metric.idx %% 1 != 0)
    stop("'metric.idx' must be integer (i.e. the index of a df column)")
  if (!(all(type.sel %in% c("mode", "median", "greater"))))
    stop("'type.sel' must be one of 'mode', 'median', 'greater'")
  if (!(all(npred.sel %in% c("min", "rnd"))))
    stop("'npred.sel' must be one of 'min', 'rnd'")


  # check the presence of NA or Inf
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(as.matrix(df))))
    stop("Inf values are not allowed in the 'df' matrix")


  # function to implement the mode
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }


  df <- df[,c(metric.idx,npred.idx)]
  df <- round(df,2)
  df$orig_idx <- 1:nrow(df)

  if(type.sel == "mode" & npred.sel == "min"){
    idx_opt_model_1 <- which(df[,metric.idx] == getmode(df[,metric.idx]))
    df_1 <- df[idx_opt_model_1,]
    idx_opt_model <- df_1[which(df_1[,npred.idx] == min(df_1[,npred.idx])),3]
    idx_opt_model <- idx_opt_model[sample(length(idx_opt_model),1)]

  }else if(type.sel == "mode" & npred.sel == "rnd"){
    idx_opt_model <- which(df[,metric.idx] == getmode(df[,metric.idx]))
    idx_opt_model <- idx_opt_model[sample(length(idx_opt_model),1)]

  }else if(type.sel == "median" & npred.sel == "min"){
    idx_opt_model_1 <- which(df[,metric.idx] == median(df[,metric.idx]))
    df_1 <- df[idx_opt_model_1,]
    idx_opt_model <- df_1[which(df_1[,npred.idx] == min(df_1[,npred.idx])),3]
    idx_opt_model <- idx_opt_model[sample(length(idx_opt_model),1)]

  }else if(type.sel == "median" & npred.sel == "rnd"){
    idx_opt_model <- which(df[,metric.idx] == median(df[,metric.idx]) )
    idx_opt_model <- idx_opt_model[sample(length(idx_opt_model),1)]

  }else if(type.sel == "greater" & npred.sel == "min"){
    if (length((which(df[,metric.idx] >= th.sel))) == 0)
      stop("no models selected")
    idx_opt_model_1 <- which(df[,metric.idx] >= th.sel)
    df_1 <- df[idx_opt_model_1,]
    idx_opt_model <- df_1[which(df_1[,npred.idx] == min(df_1[,npred.idx])),3]
    idx_opt_model <- idx_opt_model[sample(length(idx_opt_model),1)]

  }else if(type.sel == "greater" & npred.sel == "rnd"){
    if (length((which(df[,metric.idx] >= th.sel))) == 0)
      stop("no models selected")
    idx_opt_model <- which(df[,metric.idx] >= th.sel )
    idx_opt_model <- idx_opt_model[sample(length(idx_opt_model),1)]

  } else{
    stop("this type.sel does not exist")
  }


  # Plot results
  data_summary <- function(data, varname, groupnames){
    summary_func <- function(x, col){
      c(counts = length(x[[col]]))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    # data_sum <- rename(data_sum, c("counts" = varname))
    return(data_sum)
  }


  df.plot <- data_summary(df,
                          colnames(df)[1],
                          c(colnames(df)[1],
                            colnames(df)[2]))
  colnames(df.plot) <- c("Metrics", "N.predictors","Counts")

  print(ggplot(aes(x = Metrics,y = N.predictors), data = df.plot) +
          geom_point(aes(alpha=Counts, size=Counts), color="blue", alpha=0.2) +
          theme_bw() +
          scale_size(range = c(5, 15)) +
          geom_point(aes(x = df[idx_opt_model,1], y = df[idx_opt_model,2]),
                     shape=120, size=8, color="red") +
          scale_y_continuous(limits = c(min(df.plot$N.predictors),
                                        max(df.plot$N.predictors)+1)) +
          scale_x_continuous(limits = c(0.5,1),
                             breaks = c(0.5,0.55,0.6,0.65,0.7,0.75,
                                        0.8,0.85,0.9,0.95,1),
                             labels = c("0.5","0.55","0.6","0.65","0.7","0.75",
                                        "0.8","0.85","0.9","0.95","1")) +
          ggtitle("Bubble Chart")
  )

  #
  # # Print histograms
  #   print(ggplot(aes(df[,metric.idx]),data=df) +
  #           geom_histogram(aes(y=..density..),
  #                          colour="black",
  #                          fill="blue",
  #                          alpha=0.2,
  #                          breaks=seq(0.5, 1, by = 0.01))+
  #           geom_density(alpha=0.4,
  #                        fill="orange") +
  #           # geom_vline(aes(xintercept=df[90,1]),
  #           #            color="red",
  #           #            linetype="dashed", size=2)+
  #           theme_bw()+
  #           xlab(colnames(df)[metric.idx])+
  #           ylab("Counts") +
  #           ggtitle(paste0("Histogram of models performance: by ",
  #                          colnames(df)[metric.idx]))
  #   )
  #
  #
  cat("In your df, the 'optimal model' has index:",idx_opt_model, "\n")
  cat(colnames(df)[metric.idx], "=",df[idx_opt_model,1],"\n")
  cat(df[idx_opt_model,2], "predictors","\n")


  return(idx_opt_model)
}
