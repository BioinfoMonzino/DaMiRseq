## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex2(relative.path = TRUE)

## ----knitr, echo=FALSE, results="hide"-------------------------------------
library("knitr")
opts_chunk$set(
  tidy=FALSE,
  dev="png",
  fig.show="hide",
#  fig.width=4, fig.height=4.5,
  fig.width=10, fig.height=8,
  fig.pos="tbh",
  cache=TRUE,
  message=FALSE)

