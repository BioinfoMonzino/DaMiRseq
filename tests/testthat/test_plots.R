# sample dataset
data(data_norm)
data(df)
data_norm_test <- assay(data_norm[1:100,])
covar_test <- df


# check 'class' label
expect_true('class' %in% colnames(covar_test))


# expected error for DaMiR.Allplot
wrong_covar <- c("a","b")

expect_error(DaMiR.Allplot())
expect_error(DaMiR.Allplot(data_norm_test))
expect_error(DaMiR.Allplot(data_norm_test,wrong_covar))
expect_error(DaMiR.Allplot("data_norm_test",wrong_covar))
expect_error(DaMiR.Allplot(data_norm_test,covar_test,type = "foo"))

# expected error for DaMiR.MDSplot

expect_error(DaMiR.MDSplot())
expect_error(DaMiR.MDSplot(data_norm_test))
expect_error(DaMiR.MDSplot("data_norm_test"))
expect_error(DaMiR.MDSplot(data_norm_test,covar_test,type = "foo"))

# expected error for DaMiR.Clustplot
wrong_covar <- c("a","b")
data_norm_test <- t(data_norm_test)
expect_error(DaMiR.Clustplot())
expect_error(DaMiR.Clustplot(data_norm_test,wrong_covar))
expect_error(DaMiR.Clustplot("data_norm_test",wrong_covar))
expect_error(DaMiR.Clustplot(data_norm_test, covar_test,
                             type_row = "foo", type_col = "foo"))
