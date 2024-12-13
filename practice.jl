print("Hello! ")

print("how does this work")

2+2


using DataFrames
using Arrow

features_tbl = Arrow.read("data/features_20190529.feather")

first(features, 5)
features_df = DataFrame(features_tbl)
