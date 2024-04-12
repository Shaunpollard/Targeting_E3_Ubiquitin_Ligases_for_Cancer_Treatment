
# Targeting Lung Cancer 
# Abstract/Summary: 

# load packages

# List of packages
pkgs <- c("dplyr", "ComplexHeatmap", "PharmacoGx", "magicaxis",
          "tidyverse", "reshape2", "magrittr", "circlize", "plotly", 
          "gghighlight", "ggrepel", "drc", "ggpubr", "viridis")

# Load packages using lapply
lapply(pkgs, require, character.only = TRUE)



# load CRISPR knockout data of E3 ligaeses from DepMap. 
# read in data
crispr <- read.delim(file = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/metadata/CRISPR_(DepMap_Public_23Q4+Score,_Chronos)_subsetted.csv", 
                     header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE) 
crispr <- crispr[, -c(1, 4, 5, 6, 7, 8)] # get rid of unneeded columns 
# removes columns from the "crispr" data frame that contain any missing values.
crispr <- crispr[,!apply(is.na(crispr), 2, any)] 


#EGFR Mutant Lung
EGFRmutLung <- data.frame("cell_line_display_name" = c("NCIH1838", "ETCC016", "NCIH1573", "HCC827GR5", 
                                                       "NCAH3255", "HCC827", "NCAH1975", "HCC2935"), 
                          "lineage_1" = rep("Lung Mutant", 8))
crispr$lineage_1 <- ifelse(crispr$cell_line_display_name %in% EGFRmutLung$cell_line_display_name, 
                           "Lung Mutant", crispr$lineage_1)


# Count of Cancer types/linages per cell line 
rownames(crispr) <- crispr$cell_line_display_name
crispr$cell_line_display_name <- NULL
crispr_lineages_count <- crispr %>% group_by(lineage_1) %>% tally()


# Filtering the data so that we can preform a t-test on linages that have
# more than 4 samples (cell lines). I chose 4 b/c the EGRF lines only have
# 4 cell lines
crispr_lineages_count <- filter(crispr_lineages_count, n >= 4) # greater than 4 samples
#crispr_lineages_count
crispr <- subset(crispr, lineage_1 %in% crispr_lineages_count$lineage_1)


# Generating colors and breaks for lineages
lineage_colors <- viridis(length(unique(crispr$lineage_1)))
lineage_breaks <- unique(crispr$lineage_1)

# T-TEST
ttest_results <- data.frame(Lineage = lineage_breaks)

for(i in 2:ncol(crispr)) {
  genename <- colnames(crispr)[i]
  df <- crispr[, c(1, i)]
  colnames(df) <- c("Lineage", "Chronos")
  df$Chronos <- as.numeric(df$Chronos)
  any(is.na(df$Chronos))
  ttest <- compare_means(Chronos ~ Lineage,  data = df, ref.group = ".all.", method = "t.test", p.adjust.method = "BH")
  ttest <- ttest %>% dplyr::select(group2, p.adj)
  colnames(ttest) <- c("Lineage", genename)
  ttest_results <- merge(ttest_res, ttest, by = "Lineage", all = TRUE)
}
# save T-TEST results 
#write.csv(x = ttest_results, 
#          file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/T-test_results/", 
#                        Sys.Date(), 
#                        "_CRISPR_ttest_results.csv"), 
#          row.names = FALSE)


# HEATMAP

rownames(ttest_res_RNAi_achilles_demeter2) <- ttest_res_RNAi_achilles_demeter2$Lineage
ttest_res_RNAi_achilles_demeter2$Lineage <- NULL
require(circlize)


# Define row indices for "EGFR lung mutants" and "Lung"
EGFRmutLung_cancer_row_index <- which(ttest_results$Lineage == "Lung Mutant")
Lung_cancer_row_index <- which(ttest_results$Lineage == "Lung")

# Reorder rows to have "EGFR lung mutants" at the top and "Lung" second
ordered_rownames <- c(ttest_results$Lineage[EGFRmutLung_cancer_row_index], 
                      ttest_results$Lineage[Lung_cancer_row_index],
                      ttest_results$Lineage[!(ttest_results$Lineage %in% c("Lung Mutant", "Lung"))])

# Reorder the data matrix based on the new row order
ordered_matrix <- ttest_results$Lineage[ordered_rownames, ]

# Handle NA/NaN/Inf values by imputing or removing them
ordered_matrix <- na.omit(ordered_matrix)  # Remove rows with any NA/NaN/Inf values
ordered_matrix <- na.exclude(ordered_matrix)  # Alternative: Exclude NA/NaN/Inf values from calculations


# Define row indices for "EGFR lung mutants" and "Lung"
EGFRmutLung_cancer_row_index <- which(ttest_results$Lineage == "Lung Mutant")
Lung_cancer_row_index <- which(ttest_results$Lineage == "Lung")

# Reorder rows to have "EGFR lung mutants" at the top and "Lung" second
ordered_rownames <- c(ttest_results$Lineage[EGFRmutLung_cancer_row_index], 
                      ttest_results$Lineage[Lung_cancer_row_index],
                      ttest_results$Lineage[!(ttest_results$Lineage %in% c("Lung Mutant", "Lung"))])

# Define row indices for "EGFR lung mutants" and "Lung"
EGFRmutLung_cancer_row_index <- which(ttest_results$Lineage == "Lung Mutant")
Lung_cancer_row_index <- which(ttest_results$Lineage == "Lung")

# Reorder rows to have "EGFR lung mutants" at the top and "Lung" second
ordered_rownames <- c(ttest_results$Lineage[EGFRmutLung_cancer_row_index], 
                      ttest_results$Lineage[Lung_cancer_row_index],
                      ttest_results$Lineage[!(ttest_results$Lineage %in% c("Lung Mutant", "Lung"))])

# Reorder the data matrix based on the new row order
ordered_matrix <- ttest_results[ordered_rownames, ]

# Remove the lineage column
ordered_matrix$Lineage <- NULL

# Handle NA/NaN/Inf values by imputing or removing them
ordered_matrix <- na.omit(ordered_matrix)  # Remove rows with any NA/NaN/Inf values
ordered_matrix <- na.exclude(ordered_matrix)  # Alternative: Exclude NA/NaN/Inf values from calculations






# Create the heatmap using the reordered matrix
ht_RNAi_pval <- Heatmap(
  matrix = as.matrix(ordered_matrix),
  name = "Adjusted p-value",
  color_scale <- colorRamp2(breaks = c(0, 0.025, 0.05, 1), colors = c("navy", "mediumblue", "white", "white")),
  column_order = colnames(ordered_matrix),
  row_title = "Lineage",
  #row_split = rownames(ordered_matrix),
  column_title = NA,
  #clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  heatmap_legend_param = list(legend_width = unit(3, "in"),
                              at = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
                              labels = c("0", "0.01", "0.02", "0.03", "0.04", "0.05"),
                              title_position = "topcenter",
                              legend_direction = "horizontal")
)

# Print the heatmap
ht_RNAi_pval



png("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/CRISPR_heatmap.png", height = 8, width = 60, units = "in", res = 500)
draw(ht_RNAi_pval, heatmap_legend_side = "top")
#ht_RNAi_pval
dev.off()





#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


# siRNA

# load siRNA data from DepMap. 
RNAi_achilles_demeter2 <- read.table(file = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/sh:siRNA/RNAi_(Achilles,_DEMETER2)_subsetted.csv", 
                                    header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE) # read in data
RNAi_achilles_drive <- read.table(file = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/sh:siRNA/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted.csv",
                                 header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE) # read in data
RNAi_drive <- read.table(file = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/sh:siRNA/RNAi_(DRIVE,_DEMETER2)_subsetted.csv", 
                        header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE) # read in data


# get rid of columns in the data that we do not need
RNAi_achilles_demeter2 <- as.data.frame(RNAi_achilles_demeter2[, -c(1, 4, 5, 6, 7, 8)]) 
RNAi_achilles_drive <- as.data.frame(RNAi_achilles_drive[, -c(1, 4, 5, 6, 7, 8)])
RNAi_drive <- as.data.frame(RNAi_drive[, -c(1, 4, 5, 6, 7, 8)]) 

# get rid of missing data values
RNAi_achilles_demeter2_clean <- RNAi_achilles_demeter2[,!apply(is.na(RNAi_achilles_demeter2), 2, any)] # removes columns from the "RNAi_(Achilles,_DEMETER2)_subsetted" data that contain any missing values.


# Identify empty rows in the RNAi_achilles_demeter2 dataframe
#empty_rows <- unlist(sapply(RNAi_achilles_demeter2, function(x) any(x == "")))

# Replace empty strings with NA in the selected empty rows
#RNAi_achilles_demeter2[, empty_rows] <- lapply(RNAi_achilles_demeter2[, empty_rows], function(x) ifelse(x == "", NA, x))


# rows_to_delete <- numeric(0)
# 
# for (row in 1:nrow(RNAi_achilles_demeter2)) {
#   if (any(is.na(RNAi_achilles_demeter2[row,]))) {
#     rows_to_delete <- c(rows_to_delete, row)
#   }
# }
# 
# # Remove rows with NA values
# RNAi_achilles_demeter2_clean <- RNAi_achilles_demeter2[-rows_to_delete, ]
# 








# Repeat the same process for RNAi_achilles_drive and RNAi_drive data frames
#empty_rows <- sapply(RNAi_achilles_drive, function(x) any(x == ""))
#RNAi_achilles_drive[empty_rows] <- lapply(RNAi_achilles_drive[empty_rows], function(x) ifelse(x == "", NA, x))


rows_to_delete <- numeric(0)

for (row in 1:nrow(RNAi_achilles_drive)) {
  if (any(is.na(RNAi_achilles_drive[row,]))) {
    rows_to_delete <- c(rows_to_delete, row)
  }
}

# Remove rows with NA values
RNAi_achilles_drive_clean <- RNAi_achilles_drive[-rows_to_delete, ]








# Repeat the same process for RNAi_drive
#empty_rows <- sapply(RNAi_drive, function(x) any(x == ""))
#RNAi_drive[empty_rows] <- lapply(RNAi_drive[empty_rows], function(x) ifelse(x == "", NA, x))

rows_to_delete <- numeric(0)

for (row in 1:nrow(RNAi_drive)) {
  if (any(is.na(RNAi_drive[row,]))) {
    rows_to_delete <- c(rows_to_delete, row)
  }
}

# Remove rows with NA values
RNAi_drive_clean <- RNAi_drive[-rows_to_delete, ]


# EGFR Mutant Lung
# add 'Mutant Lung' as lineage for the cell lines known to be EGFRmutLung

EGFRmutLung <- data.frame("cell_line_display_name" = c("NCIH1838", "ETCC016", 
                                                       "NCIH1573", "HCC827GR5", 
                                                       "NCAH3255", "HCC827", "NCAH1975", 
                                                       "HCC2935"), "lineage_1" = rep("Lung Mutant", 8))

RNAi_achilles_demeter2_clean$lineage_1 <- ifelse(RNAi_achilles_demeter2_clean$cell_line_display_name %in% EGFRmutLung$cell_line_display_name, 
                                                 "Lung Mutant", RNAi_achilles_demeter2_clean$lineage_1)

RNAi_achilles_drive_clean$lineage_1 <- ifelse(RNAi_achilles_drive_clean$cell_line_display_name %in% EGFRmutLung$cell_line_display_name, 
                                              "Lung Mutant", RNAi_achilles_drive_clean$lineage_1)

RNAi_drive_clean$lineage_1 <- ifelse(RNAi_drive_clean$cell_line_display_name %in% EGFRmutLung$cell_line_display_name, 
                                     "Lung Mutant", RNAi_drive_clean$lineage_1)

EGFRmutLung$cell_line_display_name[EGFRmutLung$cell_line_display_name %in% RNAi_achilles_demeter2[rows_to_delete, "cell_line_display_name"]]
# "HCC827GR5" "HCC827" 
EGFRmutLung$cell_line_display_name[EGFRmutLung$cell_line_display_name %in% RNAi_achilles_drive[rows_to_delete, "cell_line_display_name"]]
# "NCIH1838" "NCIH1573" "HCC827"  
EGFRmutLung$cell_line_display_name[EGFRmutLung$cell_line_display_name %in% RNAi_drive[rows_to_delete, "cell_line_display_name"]]
# "NCIH1838" "HCC827"



EGFRmutLung$cell_line_display_name[EGFRmutLung$cell_line_display_name %in% RNAi_achilles_demeter2_clean$cell_line_display_name]
# they are removed in the clean data set.

EGFRmutLung$cell_line_display_name[EGFRmutLung$cell_line_display_name %in% RNAi_drive_clean$cell_line_display_name]
# "NCIH1838" "NCIH1573" "HCC827"  


# Tally of the Cancer types samples

rownames(RNAi_achilles_demeter2_clean) <- RNAi_achilles_demeter2_clean$cell_line_display_name
RNAi_achilles_demeter2_clean$cell_line_display_name <- NULL
RNAi_achilles_demeter2_clean_lineages_count <- RNAi_achilles_demeter2_clean %>% group_by(lineage_1) %>% tally()

rownames(RNAi_achilles_drive_clean) <- RNAi_achilles_drive_clean$cell_line_display_name
RNAi_achilles_drive_clean$cell_line_display_name <- NULL
RNAi_achilles_drive_clean_lineages_count <- RNAi_achilles_drive_clean %>% group_by(lineage_1) %>% tally()

rownames(RNAi_drive_clean) <- RNAi_drive_clean$cell_line_display_name
RNAi_drive_clean$cell_line_display_name <- NULL
RNAi_drive_clean_lineages_count <- RNAi_drive_clean %>% group_by(lineage_1) %>% tally()


# RNAi_lineages_count save csv files
#write.csv(x = RNAi_achilles_demeter2_clean_lineages_count, file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/E3_list/", Sys.Date(), "_RNAi_achilles_demeter2_lineages_count.csv"), row.names = FALSE)
#write.csv(x = RNAi_achilles_drive_clean_lineages_count, file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/E3_list/", Sys.Date(), "_RNAi_achilles_drive_lineages_count.csv"), row.names = FALSE)
#write.csv(x = RNAi_drive_clean_lineages_count, file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/E3_list/", Sys.Date(), "_RNAi_drive_lineages_count.csv"), row.names = FALSE)


# Only use Cancer types that have more than or equal to 4 samples. Update RNAi dataframes 

RNAi_achilles_demeter2_clean_lineages_count <- filter(RNAi_achilles_demeter2_clean_lineages_count, n >= 2) 
RNAi_achilles_demeter2_clean <- subset(RNAi_achilles_demeter2_clean, lineage_1 %in% RNAi_achilles_demeter2_clean_lineages_count$lineage_1)

# not sure what to do with this data
RNAi_achilles_drive_clean_lineages_count <- filter(RNAi_achilles_drive_clean_lineages_count, n >= 5) 
RNAi_achilles_drive_clean <- subset(RNAi_achilles_drive_clean, lineage_1 %in% RNAi_achilles_drive_clean_lineages_count$lineage_1)


RNAi_drive_clean_lineages_count <- filter(RNAi_drive_clean_lineages_count, n >= 3) 
RNAi_drive_clean <- subset(RNAi_drive_clean, lineage_1 %in% RNAi_drive_clean_lineages_count$lineage_1)


# Generating colors and breaks for lineages

lineage_colors <- viridis(length(unique(RNAi_achilles_demeter2_clean$lineage_1)))
lineage_breaks <- unique(RNAi_achilles_demeter2_clean$lineage_1)

lineage_colors <- viridis(length(unique(RNAi_achilles_drive_clean$lineage_1)))
lineage_breaks <- unique(RNAi_achilles_drive_clean$lineage_1)

lineage_colors <- viridis(length(unique(RNAi_drive_clean$lineage_1)))
lineage_breaks <- unique(RNAi_drive_clean$lineage_1)


# T TEST

ttest_results <- data.frame(Lineage = lineage_breaks)

for(i in 2:ncol(RNAi_achilles_demeter2_clean)) {
  genename <- colnames(RNAi_achilles_demeter2_clean)[i]
  df_RNAi_achilles_demeter2_clean <- RNAi_achilles_demeter2_clean[, c(1, i)]
  colnames(df_RNAi_achilles_demeter2_clean) <- c("Lineage", "Chronos")
  df_RNAi_achilles_demeter2_clean$Chronos <- as.numeric(df_RNAi_achilles_demeter2_clean$Chronos)
  any(is.na(df_RNAi_achilles_demeter2_clean$Chronos))
  ttest_RNAi_drive <- compare_means(Chronos ~ Lineage,  data = df_RNAi_achilles_demeter2_clean, 
                                    ref.group = ".all.", method = "t.test", p.adjust.method = "BH")
  ttest_RNAi_drive <- ttest_RNAi_drive %>% dplyr::select(group2, p.adj)
  colnames(ttest_RNAi_drive) <- c("Lineage", genename)
  ttest_results <- merge(ttest_results, ttest_RNAi_drive, 
                                by = "Lineage", all = TRUE)
}
#write.csv(x = ttest_results, 
#file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/T-test_results/", 
#              Sys.Date(), "_RNAi_achilles_demeter2_ttest_results.csv"), row.names = FALSE)


ttest_results <- data.frame(Lineage = lineage_breaks)

for(i in 2:ncol(RNAi_achilles_drive_clean)) {
  genename <- colnames(RNAi_achilles_drive_clean)[i]
  df_RNAi_achilles_drive_clean <- RNAi_achilles_drive_clean[, c(1, i)]
  colnames(df_RNAi_achilles_drive_clean) <- c("Lineage", "Chronos")
  df_RNAi_achilles_drive_clean$Chronos <- as.numeric(df_RNAi_achilles_drive_clean$Chronos)
  any(is.na(df_RNAi_achilles_drive_clean$Chronos))
  ttest_RNAi_drive <- compare_means(Chronos ~ Lineage,  data = df_RNAi_achilles_drive_clean, 
                                    ref.group = ".all.", method = "t.test", p.adjust.method = "BH")
  ttest_RNAi_drive <- ttest_RNAi_drive %>% dplyr::select(group2, p.adj)
  colnames(ttest_RNAi_drive) <- c("Lineage", genename)
  ttest_results <- merge(ttest_results, ttest_RNAi_drive, 
                                by = "Lineage", all = TRUE)
}
#write.csv(x = ttest_results, 
#file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/E3_list/", 
#  Sys.Date(), "_RNAi_achilles_drive_ttest_results.csv"), row.names = FALSE)



ttest_results <- data.frame(Lineage = lineage_breaks)

for(i in 2:ncol(RNAi_drive_clean)) {
  genename <- colnames(RNAi_drive_clean)[i]
  df_RNAi_drive <- RNAi_drive_clean[, c(1, i)]
  colnames(df_RNAi_drive) <- c("Lineage", "Chronos")
  df_RNAi_drive$Chronos <- as.numeric(df_RNAi_drive$Chronos)
  any(is.na(df_RNAi_drive$Chronos))
  ttest_RNAi_drive <- compare_means(Chronos ~ Lineage,  data = df_RNAi_drive, 
                                    ref.group = ".all.", method = "t.test", p.adjust.method = "BH")
  ttest_RNAi_drive <- ttest_RNAi_drive %>% dplyr::select(group2, p.adj)
  colnames(ttest_RNAi_drive) <- c("Lineage", genename)
  ttest_results <- merge(ttest_results, ttest_RNAi_drive, 
                                by = "Lineage", all = TRUE)
}

#write.csv(x = ttest_results, 
#file = paste0("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/T-test_results/", 
#              Sys.Date(), "_RNAi_drive_ttest_results.csv"), row.names = FALSE)





# HEAT MAP


rownames(ttest_results) <- ttest_results$Lineage
#ttest_results$Lineage <- NULL
require(circlize)


# Define row indices for "EGFR lung mutants" and "Lung"
EGFRmutLung_cancer_row_index <- which(ttest_results$Lineage == "Lung Mutant")
Lung_cancer_row_index <- which(ttest_results$Lineage == "Lung")

# Reorder rows to have "EGFR lung mutants" at the top and "Lung" second
ordered_rownames <- c(ttest_results$Lineage[EGFRmutLung_cancer_row_index], 
                      ttest_results$Lineage[Lung_cancer_row_index],
                      ttest_results$Lineage[!(ttest_results$Lineage %in% c("Lung Mutant", "Lung"))])

# Reorder the data matrix based on the new row order
ordered_matrix <- ttest_results[match(ordered_rownames, ttest_results$Lineage), ]

# Handle NA/NaN/Inf values by imputing or removing them
ordered_matrix <- na.omit(ordered_matrix)  # Remove rows with any NA/NaN/Inf values
ordered_matrix <- na.exclude(ordered_matrix)  # Alternative: Exclude NA/NaN/Inf values from calculations



# Remove the lineage column
ordered_matrix$Lineage <- NULL

# Create the heatmap using the reordered matrix
ht_RNAi_pval <- Heatmap(
  matrix = as.matrix(ordered_matrix),
  name = "Adjusted p-value",
  color_scale <- colorRamp2(breaks = c(0, 0.025, 0.05, 1), colors = c("navy", "mediumblue", "white", "white")),
  column_order = colnames(ordered_matrix),
  row_title = "Lineage",
  #row_split = rownames(ordered_matrix),
  column_title = NA,
  #clustering_distance_rows = "pearson",
  clustering_distance_columns = "euclidean",
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  heatmap_legend_param = list(legend_width = unit(3, "in"),
                              at = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
                              labels = c("0", "0.01", "0.02", "0.03", "0.04", "0.05"),
                              title_position = "topcenter",
                              legend_direction = "horizontal")
)

# Print the heatmap
ht_RNAi_pval


#png("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/RNAi_achilles_demeter2_heatmap.png", height = 8, width = 60, units = "in", res = 500)
#draw(ht_RNAi_pval, heatmap_legend_side = "top")
#ht_RNAi_pval
#dev.off()

#png("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/RNAi_achilles_drive_heatmap.png", height = 8, width = 60, units = "in", res = 500)
#draw(ht_RNAi_pval, heatmap_legend_side = "top")
#ht_RNAi_pval
#dev.off()

#png("~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/RNAi_drive_heatmap.png", height = 8, width = 60, units = "in", res = 500)
#draw(ht_RNAi_pval, heatmap_legend_side = "top")
#ht_RNAi_pval
#dev.off()




#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

# Venn Diagram

library(VennDiagram)

crispr <- c("MDM2", "PML", "HUWE1", "KMT2C", "UBR7", "TRIML2")

siRNA1 <- c("HECTD3", "AFF4", "RNF144B", "TRIM48")

siRNA2 <- c("MDM4", "TRIM7", "RC3H2", "RNF20", "RLIM", "RNF141", "BMI1", "HUWE1", 
            "RNF167", "BIRC8", "RNF170", "PCGF5", "SCAF11", "TRIM52", "TRIM13", 
            "TRIM3", "IRF2BPL", "HERC1", "PLAG1", "TRIM49C", "HECTD2", "RAG1", "RNF123")

venn.diagram(
  x = list(crispr, siRNA1, siRNA2),
  category.names = c("CRISPR", "siRNA" , "siRNA"),
  disable.logging = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  imagetype = "png",
  compression = "lzw",
  filename = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/E3_venn_diagram1.png",
  main = "E3 Overlap"
)



#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################



ttest_res <- data.frame(Lineage = lineage_breaks)
library(ggplot2)
library(ggpubr)

# Assuming "crispr" is your data frame

# for(i in 2:ncol(crispr)) {
#   genename <- colnames(crispr)[i]
#   df <- crispr[, c(1, i)]
#   colnames(df) <- c("Lineage", "Chronos")
#   df <- df[df$Lineage %in% c("Lung Mutant", "Lung"), ]  # Filter data for specified lineages
#   anova_res <- aov(Chronos ~ Lineage,  data = df)
#   p_adjusted <- summary(anova_res)$coefficients[, "Pr(>F)"] # Extract p-values from ANOVA results
#   plot <- ggplot(data = df, mapping = aes(y = reorder(Lineage, Chronos, median), x = Chronos)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(aes(color = Lineage), alpha = 0.6, height = 0.1) +
#     stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", ) +
#     labs(title = genename, subtitle = paste0("ANOVA, p-adj = ", p_adjusted[1]), x = "Chronos Dependency Score", y = "Lineage") +
#     scale_color_manual(values = lineage_colors, breaks = lineage_breaks) +
#     theme_linedraw() +
#     theme(legend.position = "none", axis.text.y = element_text(angle = 0))
#   ggsave(filename = paste0(Sys.Date(), "_", genename, "_boxplot_test.png"), path = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer", dpi = 300, device = "png", units = "in", height = 6, width = 10)
# }









# for(i in 2:ncol(crispr)) {
#   genename <- colnames(crispr)[i]
#   df <- crispr[, c(1, i)]
#   colnames(df) <- c("Lineage", "Chronos")
#   df <- df %>% filter(Lineage == "Lung Mutant" | Lineage != "Lung")
#   anova_res <- aov(Chronos ~ Lineage,  data = df)
#   p_adjusted <- summary(anova_res)$coefficients[, "Pr(>F)"] # Extract p-values from ANOVA results
#   plot <- ggplot(data = df, mapping = aes(y = reorder(Lineage, Chronos, median), x = Chronos)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(aes(color = Lineage), alpha = 0.6, height = 0.1) +
#     stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", ) +
#     labs(title = genename, subtitle = paste0("ANOVA, p-adj = ", p_adjusted[1]), x = "Chronos Dependency Score", y = "Lineage") +
#     scale_color_manual(values = lineage_colors, breaks = lineage_breaks) +
#     theme_linedraw() +
#     theme(legend.position = "none", axis.text.y = element_text(angle = 0))
#   ggsave(filename = paste0(Sys.Date(), "_", genename, "_boxplot_test.png"), path = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/ANOVA_lung_mutant_all", dpi = 300, device = "png", units = "in", height = 6, width = 10)
# }





for(i in 2:ncol(crispr)) {
  genename <- colnames(crispr)[i]
  df <- crispr[, c(1, i)]
  colnames(df) <- c("Lineage", "Chronos")
  df <- df %>% filter(Lineage == "Lung Mutant" | Lineage != "Lung")
  anova_res <- aov(Chronos ~ Lineage, data = df)
  
  # Extract p-values from ANOVA results
  p_adjusted <- summary(anova_res)[[1]]$`Pr(>F)`[1]
  
  plot <- ggplot(data = df, mapping = aes(y = reorder(Lineage, Chronos, median), x = Chronos)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Lineage), alpha = 0.6, height = 0.1) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
    labs(title = genename, subtitle = paste0("ANOVA, p-adj = ", p_adjusted), x = "Chronos Dependency Score", y = "Lineage") +
    scale_color_manual(values = lineage_colors, breaks = lineage_breaks) +
    theme_linedraw() +
    theme(legend.position = "none", axis.text.y = element_text(angle = 0))
  
  ggsave(filename = paste0(Sys.Date(), "_", genename, "_boxplot_test.png"), path = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/ANOVA_lung_mutant_all/attempt2/", dpi = 300, device = "png", units = "in", height = 6, width = 10)
}










for(i in 2:ncol(crispr)) {
  genename <- colnames(crispr)[i]
  df <- crispr[, c(1, i)]
  colnames(df) <- c("Lineage", "Chronos")
  
  # Filter the data to have two groups: Mutant EGFR Lung and all other tumors
  df <- df %>% mutate(Group = ifelse(Lineage == "Lung Mutant", "Mutant EGFR Lung", "Other Tumors"))
  
  # Perform ANOVA
  anova_res <- aov(Chronos ~ Group, data = df)
  
  # Extract p-values from ANOVA results
  p_adjusted <- summary(anova_res)[[1]]$`Pr(>F)`[1]
  
  # Plot
  plot <- ggplot(data = df, mapping = aes(y = reorder(Group, Chronos, median), x = Chronos)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Group), alpha = 0.6, height = 0.1) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
    labs(title = genename, subtitle = paste0("ANOVA, p-adj = ", p_adjusted), x = "Chronos Dependency Score", y = "Group") +
    scale_color_manual(values = c("Mutant EGFR Lung" = "blue", "Other Tumors" = "red")) +
    theme_linedraw() +
    theme(legend.position = "bottom", axis.text.y = element_text(angle = 0))
  
  ggsave(filename = paste0(Sys.Date(), "_", genename, "_boxplot_test.png"), 
         path = "~/Desktop/Weill_Cornell_Graduate/Grad_School/Thesis/Lung_Cancer/ANOVA_lung_mutant_all/attempt3/", 
         dpi = 300, device = "png", units = "in", height = 6, width = 10)
}



