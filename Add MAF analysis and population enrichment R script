library(dplyr)
library(ggplot2)

setwd("/path/to/your/data")

#------------ Data Modification for Analysis ------------#

gnomad <- read.table("./all_variant.txt", sep = "\t", header = TRUE)

cols_to_keep <- grepl("variant_id|_ac|_an", colnames(gnomad), ignore.case = TRUE) &
  !grepl("_hom|_hemi", colnames(gnomad), ignore.case = TRUE)

gnomad1 <- gnomad[, cols_to_keep]

# Modify column names
old <- c("variant_id", "joint_", "ac", "an",
         "afr", "amr", "asj", "eas", "nfe",
         "fin", "mid", "sas", "ami", "remaining")

new <- c("VCF", "", "AC", "AN",
         "AFR", "AMR", "ASJ", "EAS", "NFE",
         "FIN", "MID", "SAS", "AMI", "OTH")

colnames(gnomad1) <- sapply(colnames(gnomad1), function(x) {
  for (i in seq_along(old)) x <- gsub(old[i], new[i], x)
  x
})

gnomad1$VCF <- paste0("chr", gnomad1$VCF)

#------------ Merge with variant list ------------#

variant_list <- read.table("./variant_list.txt", sep = "\t", header = TRUE)

filtered_data <- merge(x = gnomad1,
                       y = variant_list[, c("VCF", "AC", "AN", "ACMG")],
                       by = "VCF",
                       all.y = TRUE)

#------------ Prepare AN absent columns ------------#

add_subtracted_columns <- function(data) {
  start_indices <- seq(3, 25, by = 2)
  for (i in start_indices) {
    if (i <= ncol(data)) {
      new_col_name <- paste0(colnames(data)[i], "_ab")
      data[[new_col_name]] <- data[[i]] - data[[i - 1]]
    }
  }
  return(data)
}

ft_data <- add_subtracted_columns(filtered_data)

#------------ Reorder Columns ------------#

cols <- c("VCF", "AC", "AN", "AN_ab",
          "AFR_AC", "AFR_AN", "AFR_AN_ab",
          "AMR_AC", "AMR_AN", "AMR_AN_ab",
          "ASJ_AC", "ASJ_AN", "ASJ_AN_ab",
          "EAS_AC", "EAS_AN", "EAS_AN_ab",
          "FIN_AC", "FIN_AN", "FIN_AN_ab",
          "MID_AC", "MID_AN", "MID_AN_ab",
          "NFE_AC", "NFE_AN", "NFE_AN_ab",
          "AMI_AC", "AMI_AN", "AMI_AN_ab",
          "SAS_AC", "SAS_AN", "SAS_AN_ab",
          "OTH_AC", "OTH_AN", "OTH_AN_ab")

indices <- match(cols, colnames(ft_data))
ft_data_reordered <- ft_data[, indices]

df <- ft_data_reordered[, c(-3, -6, -9, -12, -15, -18,
                            -21, -24, -27, -30, -33, -36)]

#------------ Fisher Exact Test ------------#

f <- df[, 2:25]
rownames(f) <- df$VCF

column_pairs <- list(
  ALL = c(1, 2, 3, 4),
  AFR = c(1, 2, 5, 6),
  AMR = c(1, 2, 7, 8),
  ASJ = c(1, 2, 9, 10),
  EAS = c(1, 2, 11, 12),
  FIN = c(1, 2, 13, 14),
  MID = c(1, 2, 15, 16),
  NFE = c(1, 2, 17, 18),
  AMI = c(1, 2, 19, 20),
  SAS = c(1, 2, 21, 22),
  OTH = c(1, 2, 23, 24)
)

for (group in names(column_pairs)) {
  f[[paste0("gnomad_", group, "_bonf")]] <- sapply(1:nrow(f), function(i) {
    vals <- as.numeric(f[i, column_pairs[[group]]])
    if (length(vals) == 4 && all(!is.na(vals)) && all(vals >= 0)) {
      mat <- matrix(vals, nrow = 2)
      if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
        return(tryCatch(fisher.test(mat)$p.value,
                        error = function(e) NA))
      }
    }
    NA
  }) |> p.adjust(method = "bonferroni")
}

#------------ MAF Calculation ------------#

population_pairs <- list(
  c("AC", "AN", "ALL_MAF"),
  c("AFR_AC", "AFR_AN", "AFR_MAF"),
  c("AMR_AC", "AMR_AN", "AMR_MAF"),
  c("ASJ_AC", "ASJ_AN", "ASJ_MAF"),
  c("EAS_AC", "EAS_AN", "EAS_MAF"),
  c("FIN_AC", "FIN_AN", "FIN_MAF"),
  c("MID_AC", "MID_AN", "MID_MAF"),
  c("NFE_AC", "NFE_AN", "NFE_MAF"),
  c("AMI_AC", "AMI_AN", "AMI_MAF"),
  c("SAS_AC", "SAS_AN", "SAS_MAF"),
  c("OTH_AC", "OTH_AN", "OTH_MAF")
)

filtered_data[is.na(filtered_data)] <- 0

for (pair in population_pairs) {
  filtered_data[[pair[3]]] <- filtered_data[[pair[1]]] /
    ifelse(filtered_data[[pair[2]]] == 0,
           1,
           filtered_data[[pair[2]]])
}

#------------ Long Format ------------#

maf_cols <- grep("_MAF$", names(filtered_data), value = TRUE)
bonf_cols <- grep("gnomad_.*_bonf$", names(f), value = TRUE)

maf_long <- do.call(rbind, lapply(maf_cols, function(col) {
  data.frame(VCF = filtered_data$VCF,
             Population = col,
             MAF = filtered_data[[col]])
}))

bonf_long <- do.call(rbind, lapply(bonf_cols, function(col) {
  pop <- sub("^gnomad_", "", sub("_bonf$", "_MAF", col))
  data.frame(VCF = filtered_data$VCF,
             Population = pop,
             f_pvalue = f[[col]])
}))

final_data <- merge(maf_long, bonf_long,
                    by = c("VCF", "Population"))

# Add ACMG
master <- read.table("./variant_list.txt",
                     header = TRUE,
                     sep = "\t")

acmg <- merge(master[, c("VCF", "ACMG")],
              final_data,
              by = "VCF",
              all.y = TRUE)

# Order populations
level <- c("ALL_MAF", "AFR_MAF", "AMR_MAF",
           "ASJ_MAF", "EAS_MAF", "FIN_MAF", "MID_MAF",
           "NFE_MAF", "AMI_MAF", "SAS_MAF", "OTH_MAF")

acmg$Population <- factor(acmg$Population, levels = level)
acmg$VCF <- factor(acmg$VCF, levels = unique(acmg$VCF))

#------------ Plot ------------#

png("./Figure.png",
    height = 2800,
    width = 4000,
    units = "px",
    res = 300,
    pointsize = 12)

ggplot(acmg,
       aes(x = VCF,
           y = Population,
           colour = ACMG,
           size = MAF)) +
  geom_point() +
  geom_point(data = subset(acmg, f_pvalue < 0.05),
             shape = 1,
             colour = "#D36D68",
             size = 6) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 14,
                                  face = "bold",
                                  hjust = 0)) +
  labs(title = "Population-Specific Variant Frequency Enrichment",
       y = "Population")

dev.off()
