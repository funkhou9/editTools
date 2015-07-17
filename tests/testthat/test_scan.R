library(editTools)
context("Test RNA edit search functionality")

LD_plus <- read_vcf("LD_plus_sample.vcf.gz", 2, c("DNA", "RNA"))
LD_minus <- read_vcf("LD_minus_sample.vcf.gz", 2, c("DNA", "RNA"))
liver_plus <- read_vcf("liver_plus_sample.vcf.gz", 2, c("DNA", "RNA"))
liver_minus <- read_vcf("liver_minus_sample.vcf.gz", 2, c("DNA", "RNA"))
fat_plus <- read_vcf("fat_plus_sample.vcf.gz", 2, c("DNA", "RNA"))
fat_minus <- read_vcf("fat_minus_sample.vcf.gz", 2, c("DNA", "RNA"))

test_that("Edit summary finds expected sites", {
  expect_equal(edit_summary(liver_plus)$Calls[, "POS"], c(109593,
                                                          110920,
                                                          114093,
                                                          154141))
  expect_equal(edit_summary(liver_plus)$Calls[, "CHROM"], c(10,
                                                            10,
                                                            10,
                                                            10))
  
  expect_equal(edit_summary(fat_plus)$Calls[, "POS"], c(114093,
                                                        163291))
  expect_equal(edit_summary(fat_plus)$Calls[, "CHROM"], c(10,
                                                          10))
})





