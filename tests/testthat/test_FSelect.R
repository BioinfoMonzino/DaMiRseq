## create sample dataset
samples_num <- 30
genes_num <- 1000

data_norm_test <- matrix(log2(rnbinom(genes_num*samples_num,
                                      mu=150,
                                      size=2)),
                         ncol=samples_num)
rownames(data_norm_test)<- paste0("gene", 1:genes_num)
data_norm_test[is.infinite(data_norm_test)] <- 5

covar_test <- data.frame(class=c(rep("A",samples_num/2),
                                 rep("B",samples_num/2)),
                         row.names = paste0("V", 1:samples_num))

data_norm_test <- as.data.frame(t(data_norm_test))
rownames(data_norm_test) <- rownames(covar_test)

####################################################################
## test
# check 'class' label
expect_true('class' %in% colnames(covar_test))

#expected error
expect_error(DaMiR.FSelect())
expect_error(DaMiR.FSelect("character","character"))
expect_error(DaMiR.FSelect(data_norm_test,
                           covar_test,
                           th.corr = "character"))
expect_error(DaMiR.FSelect(data_norm_test,
                           covar_test,
                           type = "foo"))
expect_error(DaMiR.FSelect(data_norm_test,
                           covar_test,
                           th.VIP =  "character"))

#launch script
testOut <- DaMiR.FSelect(data_norm_test,
                         covar_test,
                         th.corr=0.01,
                         th.VIP=0.01)

# check results
expect_true(dim(testOut$data)[1] == dim(data_norm_test)[1])
expect_true(dim(testOut$data)[2] <= dim(data_norm_test)[2])

# too stringent filtering
expect_null(testOut <- DaMiR.FSelect(data_norm_test,
                                     covar_test,
                                     th.corr=0.01,
                                     th.VIP = 10))
