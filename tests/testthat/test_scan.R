library(ErnstEditTools)
context("Test the scanning for edits - one step at a time")

LD_plus <- read_vcf("LD_plus_sample.vcf.gz", 2, c("DNA", "RNA"))
LD_minus <- read_vcf("LD_minus_sample.vcf.gz", 2, c("DNA", "RNA"))
liver_plus <- read_vcf("liver_plus_sample.vcf.gz", 2, c("DNA", "RNA"))
liver_minus <- read_vcf("liver_minus_sample.vcf.gz", 2, c("DNA", "RNA"))
fat_plus <- read_vcf("fat_plus_sample.vcf.gz", 2, c("DNA", "RNA"))
fat_minus <- read_vcf("fat_minus_sample.vcf.gz", 2, c("DNA", "RNA"))




liver_plus_edits <- edit_summary(liver_plus)

fat_plus_edits <- edit_summary(fat_plus)
fat_minus_edits <- edit_summary(fat_minus)


edit_summary(fat_plus)
