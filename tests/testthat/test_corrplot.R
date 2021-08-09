#sample dataset
data(data_norm)
data(sv)
df <- colData(data_norm)
#DaMiR.corrplot(sv, df, type="pearson", sig.level=0.0001)

#expected error
#expect_error(DaMiR.corrplot())
#expect_error(DaMiR.corrplot(sv))
#expect_error(DaMiR.corrplot("sv"))
#expect_error(DaMiR.corrplot(sv, df, type="foo"))
#expect_error(DaMiR.corrplot(sv, sig.level="character"))

