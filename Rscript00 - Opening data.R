# =========================================================================
# OPENING DATASET IN R
#
# This section loads the pooled respondent-by-debate data used in the
# subsequent analysis. (see main text for pooled sample composition in
# Chapter 4.2)
#
# The loaded object is stored as data_pooled and serves as the main input
# dataset for the analyses that follow.
# =========================================================================

library(dplyr)
library(forcats)

# Download data_pooled.rds and insert the local file path below.

data_pooled <- readRDS("C:/.............../data_pooled.rds")

