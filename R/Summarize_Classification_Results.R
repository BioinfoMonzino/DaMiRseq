
#' Plot Classification results
#'
#' @param Classification_res a list returned from DaMiR.EnsembleLearning2cl
#'  or DaMiR.EnsembleLearningNcl
#'
#' @return the requested plot
#'
DaMiR.EnsembleLearningPlot <- function(Classification_res,what= c("accuracy","MCC","Specif",
                                                                  "Sensit")){

  if (missing(what)){
    what <- names(Classification_res)
  }
  if ("accuracy" %in% what){
    acc_dotplot <- melt(as.data.frame(Classification_res$accuracy),
                        measure.vars = colnames(Classification_res$accuracy))

    colnames(acc_dotplot) <- c("Classifiers","Accuracy")
    print(ggplot(acc_dotplot, aes(x=Classifiers,y=Accuracy)) +
            geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
            geom_dotplot(binaxis='y',
                         stackdir='center',
                         stackratio=1.5,
                         dotsize=0.2,
                         binwidth = 0.5) +
            stat_summary(fun.data=mean_sdl,
                         fun.args = list(mult=1),
                         geom="pointrange",
                         color="white") +
            coord_cartesian(ylim=c(min(Classification_res$accuracy)-5,100))

    )
  }

  if ("MCC" %in% what){
    mcc_dotplot <- melt(as.data.frame(Classification_res$MCC),
                        measure.vars = colnames(Classification_res$MCC))

    colnames(mcc_dotplot) <- c("Classifiers","MCC")
    print(ggplot(mcc_dotplot, aes(x=Classifiers,y=MCC)) +
            geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
            geom_dotplot(binaxis='y',
                         stackdir='center',
                         stackratio=1.5,
                         dotsize=0.002,
                         binwidth = 0.5) +
            stat_summary(fun.data=mean_sdl,
                         fun.args = list(mult=1),
                         geom="pointrange",
                         color="white") +
            coord_cartesian(ylim=c(min(Classification_res$MCC)-0.05,1))
    )
  }
  if ("Specif" %in% what){

    spe_dotplot <- melt(as.data.frame(Classification_res$Specif),
                        measure.vars = colnames(Classification_res$Specif))

    colnames(spe_dotplot) <- c("Classifiers","Specificity")
    print(ggplot(spe_dotplot, aes(x=Classifiers,y=Specificity)) +
            #ylim(0,1) +
            geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
            geom_dotplot(binaxis='y',
                         stackdir='center',
                         stackratio=1.5,
                         dotsize=0.002,
                         binwidth = 0.5) +
            stat_summary(fun.data=mean_sdl,
                         fun.args = list(mult=1),
                         geom="pointrange",
                         color="white") +
            coord_cartesian(ylim=c(min(Classification_res$Specif)-0.05,1))
    )
  }
  if ("Sensit" %in% what){
    sen_dotplot <- melt(as.data.frame(Classification_res$Sensit),
                        measure.vars = colnames(Classification_res$Sensit))

    colnames(sen_dotplot) <- c("Classifiers","Sensitivity")
    print(ggplot(sen_dotplot, aes(x=Classifiers,y=Sensitivity)) +
            #ylim(0,1) +
            geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
            geom_dotplot(binaxis='y',
                         stackdir='center',
                         stackratio=1.5,
                         dotsize=0.002,
                         binwidth = 0.5) +
            stat_summary(fun.data=mean_sdl,
                         fun.args = list(mult=1),
                         geom="pointrange",
                         color="white") +
            coord_cartesian(ylim=c(min(Classification_res$Sensit)-0.05,1))
    )

  }
  # TODO use facets to plot all in one page?


}


#' Print Classification results
#'
#' @param Classification_res a list returned from DaMiR.EnsembleLearning2cl
#'  or DaMiR.EnsembleLearningNcl
#'
#' @return the requested plot
#'
DaMiR.EnsembleLearningPrint <- function(Classification_res,
                                        what= c("accuracy","MCC","Specif",
                                                "Sensit")){
  if (missing(what)){
    what <- names(Classification_res)
  }
  if ("accuracy" %in% what){
    cat("Accuracy [%]:",
        "\n",
        colnames(Classification_res$accuracy),
        "\n",
        "Mean:",round(colMeans(Classification_res$accuracy),2),"\n","St.Dev.",
        round(colSds(Classification_res$accuracy),digits = 2),"\n")}

  if ("MCC" %in% what){
    cat("MCC score:",
        "\n",
        colnames(Classification_res$MCC),
        "\n",
        "Mean:",round(colMeans(Classification_res$MCC),2),"\n","St.Dev.",
        round(colSds(Classification_res$MCC),digits = 2),"\n")}


  if ("Specif" %in% what){
    cat("Specificity:",
        "\n",
        colnames(Classification_res$Specif),
        "\n",
        "Mean:",round(colMeans(Classification_res$Specif),2),"\n","St.Dev.",
        round(colSds(Classification_res$Specif),digits = 2),"\n")}

  if ("Sensit" %in% what){
    cat("Sensitivity:",
        "\n",
        colnames(Classification_res$Sensit),
        "\n",
        "Mean:",round(colMeans(Classification_res$Sensit),2),"\n","St.Dev.",
        round(colSds(Classification_res$Sensit),digits = 2),"\n")}


  #  cat("PPV:",
  #      "\n",
  #      colnames(PPV.Class),
  #      "\n",
  #      "Mean:",round(colMeans(PPV.Class),2),"\n","St.Dev.",
  #      round(colSds(PPV.Class),digits = 2),"\n")
  #  cat("NPV:",
  #      "\n",
  #      colnames(NPV.Class),
  #      "\n",
  #      "Mean:",round(colMeans(NPV.Class),2),"\n","St.Dev.",
  #      round(colSds(NPV.Class),digits = 2),"\n")

}



