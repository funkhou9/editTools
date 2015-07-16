# Need to troubleshoot - why doesn't this run properly using 'check'?


# library(ErnstEditTools)
# context("Test vcf attributes and methods")
# 
# test_that("vcf class has working getters and fields", {
#   expect_is(read_vcf("fat_minus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               snps(), 
#             "data.frame")
#   expect_is(read_vcf("fat_plus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               snps(), 
#             "data.frame")
#   expect_is(read_vcf("liver_minus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               snps(), 
#             "data.frame")
#   expect_is(read_vcf("liver_plus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               snps(), 
#             "data.frame")
#   expect_is(read_vcf("LD_minus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               snps(), 
#             "data.frame")
#   expect_is(read_vcf("LD_plus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               snps(), 
#             "data.frame")
#   expect_is(read_vcf("fat_minus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               samples(), 
#             "list")
#   expect_is(read_vcf("fat_plus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               samples(), 
#             "list")
#   expect_is(read_vcf("liver_minus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               samples(), 
#             "list")
#   expect_is(read_vcf("liver_plus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               samples(), 
#             "list")
#   expect_is(read_vcf("LD_minus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               samples(), 
#             "list")
#   expect_is(read_vcf("LD_plus_sample.vcf.gz", 2, c("DNA", "RNA")) %>%
#               samples(), 
#             "list")
# })