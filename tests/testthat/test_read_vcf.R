library(ErnstEditTools)
context("Test the ability to read in vcf files with read_vcf()")

test_that("read_vcf() returns vcf class", {
  expect_is(read_vcf("fat_minus_sample.vcf.gz", 2, c("DNA", "RNA")),
            "vcf")
  expect_is(read_vcf("fat_plus_sample.vcf.gz", 2, c("DNA", "RNA")),
            "vcf")
  expect_is(read_vcf("liver_minus_sample.vcf.gz", 2, c("DNA", "RNA")),
            "vcf")
  expect_is(read_vcf("liver_plus_sample.vcf.gz", 2, c("DNA", "RNA")),
            "vcf")
  expect_is(read_vcf("LD_minus_sample.vcf.gz", 2, c("DNA", "RNA")),
            "vcf")
  expect_is(read_vcf("LD_plus_sample.vcf.gz", 2, c("DNA", "RNA")),
            "vcf")
})