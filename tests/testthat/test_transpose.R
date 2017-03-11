#sample dataset
data(data_norm)
testdata<-assay(data_norm)

#expected error
expect_error(DaMiR.transpose())
expect_error(DaMiR.transpose("data_norm"))

#launch function
trans_testdata <- DaMiR.transpose(testdata)

# check results
expect_true(dim(testdata)[2] == dim(trans_testdata)[1])
expect_true(dim(testdata)[1] == dim(trans_testdata)[2])
expect_is(trans_testdata, "data.frame")
