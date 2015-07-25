library(editTools)
context("Test RNA edit scan functionality")

# liver_edits <- edit_scan("liver_plus_sample.vcf.gz", "liver_minus_sample.vcf.gz")
# 
# test_that("Edit summary finds expected sites", {
#   expect_equal(edit_summary(liver_plus)$Calls[, "POS"], c(109593,
#                                                           110920,
#                                                           114093,
#                                                           154141))
#   expect_equal(edit_summary(liver_plus)$Calls[, "CHROM"], c(10,
#                                                             10,
#                                                             10,
#                                                             10))
#   
#   expect_equal(edit_summary(fat_plus)$Calls[, "POS"], c(114093,
#                                                         163291))
#   expect_equal(edit_summary(fat_plus)$Calls[, "CHROM"], c(10,
#                                                           10))
# })
# 
# t <- read_vcf("liver_plus_sample.vcf.gz", c("DNA", "RNA"))




system.time({t <- 
  editTools::read_vcf("/Users/sfunkhouser/Programs/Cpp/editTools/1502_DNA_LD_plus_all.vcf") %>%
    editTools::edit_summary(10, TRUE, 10, 95, 5, TRUE)})