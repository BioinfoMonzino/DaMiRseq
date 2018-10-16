# sample dataset
samples_num <- 11
genes_num <- 100

dataTestSV <- matrix(rnbinom(genes_num*samples_num,mu=150,size=2),ncol=samples_num)
rownames(dataTestSV)<- paste0("gene",1:genes_num)
colnames(dataTestSV) = c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11")
dataTestSV[is.infinite(dataTestSV)] <- 5

covTestSV <- data.frame(a=c(rep("aa",3),rep("bb",4),rep("cc",4)),
                        b=c(1,2,3,4,5,6,7,8,9,0,0),
                        class=c(rep("A",5),rep("B",6)),
                        row.names = c(
                          "V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11"))
SEtest <- DaMiR.makeSE(dataTestSV,covTestSV)

#expected error
expect_error(DaMiR.normalization())
expect_error(DaMiR.normalization(data = "SEtest"))
expect_error(DaMiR.normalization(SEtest, minCounts="character"))
expect_error(DaMiR.normalization(SEtest, fSample="character"))
expect_error(DaMiR.normalization(SEtest, th.cv="character"))

# check object type
expect_true("class" %in% colnames(covTestSV))

# launch function
SEtest_norm <- DaMiR.normalization(SEtest)

# check results
expect_true(dim(SEtest_norm)[1] <= dim(SEtest)[1])
expect_true(dim(SEtest_norm)[2] == dim(SEtest)[2])
expect_true(all(colnames(SEtest) == colnames(SEtest_norm)))
expect_is(SEtest_norm, "SummarizedExperiment")
expect_true("class" %in% colnames(SEtest_norm@colData))
