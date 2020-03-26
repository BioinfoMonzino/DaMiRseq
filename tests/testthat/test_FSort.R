# create sample dataset
samples_num <- 30
genes_num <- 50

data_norm_test_1 <- matrix(log2(rnbinom(genes_num*(samples_num/2),
                                        mu=1500,
                                        size=2)),
                           ncol=samples_num/2)

data_norm_test_2 <- matrix(log2(rnbinom(genes_num*(samples_num/2),
                                        mu=20,
                                        size=2)),
                           ncol=samples_num/2)


data_norm_test <- cbind(data_norm_test_1,data_norm_test_2)

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
expect_error(DaMiR.FSort())
expect_error(DaMiR.FSort("character","character"))
expect_error(DaMiR.FSort(data_norm_test,
                         covar_test,
                           fSample = "character"))

# less stringent argument
expect_error(DaMiR.FSort(data_norm_test,
                         covar_test,
                         fSample = 0))
#launch script
testOut <- DaMiR.FSort(data_norm_test,
                       covar_test)

# check results
expect_true(dim(testOut)[2] == 2)
expect_true(all(testOut$RReliefF<=1))
