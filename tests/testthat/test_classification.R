## create sample dataset
samples_num <- 30
genes_num <- 10

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
expect_error(DaMiR.EnsembleLearning())
expect_error(DaMiR.EnsembleLearning("character","character"))
expect_error(DaMiR.EnsembleLearning("character","character","character"))
expect_error(DaMiR.EnsembleLearning(data_norm_test,
                                    covar_test,
                                    fSample.tr="character"))
expect_error(DaMiR.EnsembleLearning(data_norm_test,
                                    covar_test,
                                    fSample.tr.w="character"))
expect_error(DaMiR.EnsembleLearning(data_norm_test,
                                    covar_test,
                                    iter="character"))

#launch script
set.seed(123)
iter_set <-1
classes <-covar_test$class
testOut<-DaMiR.EnsembleLearning(data_norm_test, classes, iter=iter_set)
# check results
expect_true(dim(testOut$accuracy)[1] == iter_set)
expect_true(dim(testOut$accuracy)[2] == 7)
expect_true(all(testOut$accuracy<=100))
expect_true(all(testOut$accuracy>=0))


