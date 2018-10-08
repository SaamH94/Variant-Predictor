#function for predicting whether a single nucleotide varaint in RNA seq data is SNP or RNA editing site
#load vcf file with bedr package
#separate columns containing genotype (GT) and allele depth (AD) for each allele
#calculate ratio of ref allele to alt allele
#annotate each variant based on that and genotype as editing site or SNP and discard uninformative reads


variant_predictor <- function(vcf_file, heterozygous_variant_filename_csv, uninformative_variants_filename_csv, homozygous_unknown_variants_filename_csv){
  library(bedr)
  library(stringr)
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  vcf_file <- read.vcf(vcf_file)
  vcf_data <- vcf_file$vcf
  colnames(vcf_data)[10] <- "SNP_Editing_site_info"
  dbsnps <- which(startsWith(as.character(vcf_data$ID), "rs"))
  dbsnp_data <- vcf_data[dbsnps,]
  
  
  vcf_data_split <- str_split_fixed(vcf_data$SNP_Editing_site_info, ":", 7)
  vcf_data_final <- data.frame(vcf_data, vcf_data_split)
  vcf_data_split_2 <- str_split_fixed(vcf_data_final$X2, ",", 2)
  vcf_data_final_2 <- data.frame(vcf_data$CHROM, vcf_data$POS, vcf_data$ID, vcf_data$REF, vcf_data$ALT, vcf_data$QUAL, vcf_data$FILTER, vcf_data$INFO, vcf_data_split, vcf_data_split_2)
  
  
  
  homozygous_var <- which(vcf_data_final_2$X1 == "1/1" & vcf_data_final_2$X2.1 !=0)
  homozygous_sites <- vcf_data_final_2[homozygous_var,]
  homozygous_sites$ref_alt_ratio <- as.numeric.factor(homozygous_sites$X1.1)/as.numeric.factor(homozygous_sites$X2.1)
  homozygous_sites$Variant_type = NA
  
  homozygous_sites$Variant_type[which(homozygous_sites$ref_alt_ratio>0)] <- "Editing site"
  homozygous_sites$Variant_type[which(homozygous_sites$ref_alt_ratio == 0)] <- "SNP"
  
  
  
  dbsnps_homozygous <- which(startsWith(as.character(homozygous_sites$vcf_data.ID), "rs"))
  homozygous_SNPs <- homozygous_sites[dbsnps_homozygous,]
  homozygous_non_dbsnp_pos <- which(homozygous_sites$vcf_data.ID == ".")
  homozygous_unknown_Variants <- homozygous_sites[homozygous_non_dbsnp_pos,]
  
  
  heterozygous_var <-  which(vcf_data_final_2$X1 == "0/1" & vcf_data_final_2$X2.1 != 0)
  heterozygous_sites <- vcf_data_final_2[heterozygous_var,]
  heterozygous_sites$ref_alt_ratio <- as.numeric.factor(heterozygous_sites$X1.1)/as.numeric.factor(heterozygous_sites$X2.1)
  heterozygous_sites$variant_type <- NA
  
  heterozygous_sites$X1.1 <- as.numeric.factor(heterozygous_sites$X1.1)
  heterozygous_sites$X2.1 <- as.numeric.factor(heterozygous_sites$X2.1)
  
  
  
  heterozygous_sites$variant_type[which(heterozygous_sites$ref_alt_ratio == 1)] <- "SNP"
  heterozygous_sites$variant_type[which(heterozygous_sites$ref_alt_ratio != 1)] <- "Editing site"
  heterozygous_sites$variant_type[which(heterozygous_sites$X1.1 == heterozygous_sites$X2.1)] <- "SNP"
  
  heterozygous_Uninformative_reads <- which(vcf_data_final_2$X1 == "0/1" & vcf_data_final_2$X2.1 == 0)
  heterozygous_Uninformative_reads <- vcf_data_final_2[heterozygous_Uninformative_reads,] 
  heterozygous_Uninformative_reads$ref_alt_ratio <- NA
  heterozygous_Uninformative_reads$variant_type <- "Uninformative read"
  
  homozygous_Uninformative_reads_variants <- which(vcf_data_final_2$X1 == "1/1" & vcf_data_final_2$X2.1 == 0)
  homozygous_Uninformative_reads <- vcf_data_final_2[homozygous_Uninformative_reads_variants,]
  homozygous_Uninformative_reads$ref_alt_ratio <- NA
  homozygous_Uninformative_reads$variant_type <- "Uninformative read"
  
  
  names(homozygous_sites)[19] <- "Variant_type"
  names(homozygous_Uninformative_reads)[19] <- "Variant_type"
  names(heterozygous_sites)[19] <- "Variant_type"
  names(heterozygous_Uninformative_reads)[19] <- "Variant_type"
  
  total_variants <- heterozygous_sites
  total_variants_uninformative <- rbind(heterozygous_Uninformative_reads, homozygous_Uninformative_reads)
  
  write.csv(total_variants, file = "heterozygous_variant_filename_csv")
  write.csv(total_variants_uninformative, file = "uninformative_variants_filename_csv")
  write.csv(homozygous_unknown_Variants, file = "homozygous_unknown_variants_filename_csv")
  
  
  
  
  
}



