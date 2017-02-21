# sample dataset
data(data_norm)

#expected error
expect_error(DaMiR.sampleFilt())
expect_error(DaMiR.sampleFilt("data_norm"))
expect_error(DaMiR.sampleFilt(data_norm, th.corr="character"))
expect_error(DaMiR.sampleFilt(data_norm, type="foo"))

# check object type
expect_true("class" %in% colnames(colData(data_norm)))

# launch function
SEtest_filt <- DaMiR.sampleFilt(data_norm,th.corr = 0)

# too stringent
expect_error(DaMiR.sampleFilt(data_norm, th.corr=1))

# check results
expect_true(dim(data_norm)[2] >= dim(SEtest_filt)[2])
expect_true(dim(data_norm)[1] == dim(SEtest_filt)[1])
expect_is(SEtest_filt, "SummarizedExperiment")
expect_true("class" %in% colnames(SEtest_filt@colData))
