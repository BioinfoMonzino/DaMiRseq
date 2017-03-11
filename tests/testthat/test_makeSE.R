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


# check object type

expect_true("class" %in% colnames(covTestSV))
expect_true(identical(colnames(dataTestSV), rownames(covTestSV)))

#expected error
expect_error(DaMiR.makeSE())
expect_error(DaMiR.makeSE(dataTestSV="data.frame"))
expect_error(DaMiR.makeSE(covTestSV="matrix"))

# launch function
SEtest <- DaMiR.makeSE(dataTestSV,covTestSV)
# check results
expect_true(all(rownames(SEtest) == rownames(dataTestSV)))
expect_true(all(colnames(SEtest) == colnames(dataTestSV)))
expect_is(SEtest, "SummarizedExperiment")
expect_true("class" %in% colnames(SEtest@colData))

