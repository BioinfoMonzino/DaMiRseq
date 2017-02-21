# sample dataset
data(data_norm)

# check object type
expect_true("class" %in% colnames(colData(data_norm)))


#expected error
expect_error(DaMiR.SV())
expect_error(DaMiR.SV(data_norm, method = "foo"))
expect_error(DaMiR.SV(data_norm, th.fve = "character"))

# launch function
sv <- DaMiR.SV(data_norm)

# check results
expect_true(dim(sv)[1] == dim(data_norm)[2])
expect_true(dim(sv)[2] <= dim(sv)[1])
expect_is(sv, "matrix")

