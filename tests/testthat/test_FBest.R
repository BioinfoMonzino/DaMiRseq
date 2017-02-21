samples_num <- 30
genes_num <- 15

data_norm_test <- matrix(log2(rnbinom(genes_num*samples_num,
                                      mu=150,
                                      size=2)),
                         ncol=samples_num)
rownames(data_norm_test)<- paste0("gene", 1:genes_num)

data_norm_test[is.infinite(data_norm_test)] <- 5

rTest <- data.frame(RReliefF=order(seq_len(genes_num),
                                   decreasing=TRUE)/genes_num,
                    scaled.RReliefF=scale(order(seq_len(genes_num),
                                                 decreasing=TRUE)/genes_num),
                    row.names=paste0("gene", 1:genes_num))

data_norm_test <- as.data.frame(t(data_norm_test))
############################################################
## test

#expected error
expect_error(DaMiR.FBest())
expect_error(DaMiR.FBest("character"))
expect_error(DaMiR.FBest(data_norm_test,"character"))
expect_error(DaMiR.FBest(data_norm_test,autoselect = "foo"))
expect_error(DaMiR.FBest(data_norm_test,n.pred = 60))
expect_error(DaMiR.FBest(data_norm_test,n.pred = "character"))
expect_error(DaMiR.FBest(data_norm_test,th.zscore = "character"))
expect_error(DaMiR.FBest(data_norm_test,autoselect = "yes", n.pred = 5))
expect_error(DaMiR.FBest(data_norm_test,autoselect = "no", th.zscore = 1))
expect_error(DaMiR.FBest(data_norm_test,th.zscore = "character"))

#launch script #1
n.pred_set <- 10
testOut1 <- DaMiR.FBest(data_norm_test,rTest,n.pred = n.pred_set)
# check results
expect_true(dim(testOut1$data)[2] == n.pred_set)
expect_true(dim(testOut1$data)[1] == dim(data_norm_test)[1])
expect_true(length(testOut1$predictors) == n.pred_set)

#launch script #2
th.zscore_set <-1
testOut2 <- DaMiR.FBest(data_norm_test,rTest,autoselect = "yes",th.zscore = th.zscore_set)
# check results
expect_true(dim(testOut2$data)[2] == length(which(rTest$scaled.RReliefF > th.zscore_set)))




