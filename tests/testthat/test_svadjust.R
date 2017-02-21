#sample dataset
data(data_norm)
data(sv)

#expected error
expect_error(DaMiR.SVadjust())
expect_error(DaMiR.SVadjust(data_norm))
expect_error(DaMiR.SVadjust("data_norm"))
expect_error(DaMiR.SVadjust(data_norm, sv="character"))
expect_error(DaMiR.SVadjust(data_norm, sv=as.data.frame(sv)))
expect_error(DaMiR.SVadjust(data_norm, n.sv=dim(data_norm)[2]+1))

# check object type
expect_true("class" %in% colnames(colData(data_norm)))

#launch function
data_adj<-DaMiR.SVadjust(data_norm, sv, 4)

# check results
expect_true(dim(data_norm)[2] == dim(data_adj)[2])
expect_true(dim(data_norm)[1] == dim(data_adj)[1])
expect_is(data_adj, "SummarizedExperiment")
expect_true("class" %in% colnames(data_adj@colData))
