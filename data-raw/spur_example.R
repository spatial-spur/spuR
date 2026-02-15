# Create bundled example dataset from Chetty et al. (2014) commuting zone data.
# This script reads the full fixture, selects 50 rows and key columns,
# and saves it as data/spur_example.rda.

df <- utils::read.csv(
  "tests/testthat/fixtures/chetty_data_1.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

spur_example <- df[seq_len(50), c("cz", "czname", "state", "lat", "lon",
                                   "am", "rm", "gini", "fracblack")]
rownames(spur_example) <- NULL

usethis::use_data(spur_example, overwrite = TRUE)
