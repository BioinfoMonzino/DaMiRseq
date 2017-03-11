## create sample dataset
samples_num <- 30
genes_num <- 100

data_norm_test <- matrix(log2(rnbinom(genes_num*samples_num,
                                      mu=150,
                                      size=2)),
                         ncol=samples_num)
rownames(data_norm_test)<- paste0("gene", 1:genes_num)
data_norm_test[is.infinite(data_norm_test)] <- 5
data_norm_test <- as.data.frame(t(data_norm_test))

####################################################################
## test

#expected error
expect_error(DaMiR.FReduct())
expect_error(DaMiR.FReduct("character"))
expect_error(DaMiR.FReduct(data_norm_test,
                           th.corr = "character"))
expect_error(DaMiR.FReduct(data_norm_test,
                           type = "foo"))

#launch script
testOut <- DaMiR.FReduct(data_norm_test,
                          th.corr=1)
# check results
expect_true(dim(testOut)[1] == dim(data_norm_test)[1])
expect_true(dim(testOut)[2] <= dim(data_norm_test)[2])

# too stringent filtering
expect_equal(dim(DaMiR.FReduct(data_norm_test, th.corr=1))[2],0)

