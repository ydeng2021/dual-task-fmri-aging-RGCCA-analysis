# RGCCA analysis pipeline 
# This pipeline uses motor network input files
# You can replace the input files with cognitive brain FC for cognitive network analysis
# All brain and behavior data are provided in RGCCA folder
# Author: Yan Deng
#
# Biological Psychology Lab, 
# Department of Psychology, 
# School of Medicine and Health Sciences, 
# Carl von Ossietzky Universität Oldenburg, 
# Oldenburg, Germany
# Version 3.0 updated: 2026-02-16
#
# Notes:
# - Edit only the "USER SETTINGS" section when moving machines/folders.

# Load the packages
# set the working directory 
# =========================
# Libraries
# =========================
library(dplyr)
# this library is used in data cleaning, exploration, and preparation workflows, 
# especially with the pipe operator %>% 

# Install mixOmics via Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install mixOmics
BiocManager::install("mixOmics")
# Load library mixOmics
library(mixOmics) # Run wrapper.rgcca(), Loaded mixOmics, version: 6.30.0

library(ggplot2) 
library(GGally)
library(tibble)
library(qgraph)
library(tidyr)
library(stringr)

# USER SETTINGS
# =========================
project_dir <- "Your Project Directory"
setwd(project_dir)

# Input files (edit if you do cognitive task network analysis)
fc_dual_file    <- "Old_Young_SPM_motorNetwork_Go_dual_90connections.csv"
fc_single_file  <- "Old_Young_SPM_motorNetwork_motor_single_558connections.csv"
behavior_file   <- "Merged_motor_cog_summary_data_imputed.csv"
subjects_file   <- "subjects_summary_64_young_old_clean_ID_2Groups.csv"

# Step 1: Prepare data blocks
block1 <- read.csv(fc_dual_file)      # Dual-task FC
block2 <- read.csv(fc_single_file)   # Single-task FC

# load the behavior data
block3 <- readRDS("data/block3_motor_beh_clean_20250418.rds") # Single and dual behavior 
block4 <- readRDS("data/block4_motor_beh_clean_20250418.rds") # Dual behavior only
group_vector <- readRDS("data/group_vector_clean_20250418.rds")

sum(is.na(block1))
sum(is.na(block2))
# subject 08, 14, 28, 29 excluded from brain data analysis

# Step 2: Combine into a list, 
# Behavioral measurement includes only the Dual task, 
blocks <- list(
  "DualTaskFC" = block1,
  "SingleMotorFC" = block2,
  "DualBehavior" = block4
)
# always make block names character strings, not symbols.

# Step 3: Create the design matrix
# Fully connected design (off-diagonal = 1, diagonal = 0)
# connect every block with the others (no self-links)
design <- matrix(
  c(0, 1, 1,
    1, 0, 1,
    1, 1, 0),
  nrow = 3, byrow = TRUE)

colnames(design) <- rownames(design) <- names(blocks)
# Make sure they are character strings, matching the block names.

# print and inspect the design:
design
#              DualTaskFC SingleMotorFC DualBehavior
#DualTaskFC             0             1            1
#SingleMotorFC          1             0            1
#DualBehavior           1             1            0

#################################################
#___________________ Step 4: Run wrapper.rgcca()
# Set number of components to extract per block
ncomp <- 2          # ncomp: number of components to extract per block
scheme <-"horst"    # other options:  "factorial", "centroid" (Default: "horst").
set.seed(123)       # set seed

# Run RGCCA (wrapper.rgcca)
res_rgcca_dual <- wrapper.rgcca(
  X      = blocks,
  design = design,
  ncomp  = ncomp,
  scheme = scheme,
  scale  = TRUE,
  all.outputs = TRUE
)

# Check the structure of res_rgcca_dual:
str(res_rgcca_dual, max.level = 2)
# a summary of what's inside (scores, loadings, correlations, tau values, etc.).

# 2. Component scores (aka latent variables / variates)
# # Component scores (variates) per block
res_rgcca_dual$variates
# These are the projections of each subject onto the RGCCA components.

# Read scores for difference block
head(res_rgcca_dual$variates$DualBehavior)  # Scores for the Behavior block

# 3. Loadings (aka weight vectors or a coefficients)
# tell how much each variable contributes to the component.
res_rgcca_dual$loadings$DualBehavior

# 4. top 10 features by absolute loading on component 1 for each block

# 5. Inner AVE (average variance explained between blocks)
# This reflects how well the blocks relate to each other.
res_rgcca_dual$AVE$AVE_inner

# Look at the first few component scores
head(res_rgcca_dual$variates$DualTaskFC)
head(res_rgcca_dual$variates$SingleMotorFC)

# 7 Plot individuals 
plotIndiv(res_rgcca, legend = TRUE)
# If you have a group variable (e.g., AgeGroup):
plotIndiv(res_rgcca_dual, group = group_vector, legend = TRUE)


# Plot Component Scores for All Three Blocks,
#  1. Combine scores from each block
# Extract component scores
scores_dual_fc <- res_rgcca_dual$variates$DualTaskFC
scores_motor_fc <- res_rgcca_dual$variates$SingleMotorFC
scores_behavior <- res_rgcca_dual$variates$DualBehavior

# Create data frames with block labels
df_dual_fc <- data.frame(Comp1 = scores_dual_fc[, 1],
                         Comp2 = scores_dual_fc[, 2],
                         Group = group_vector,
                         Block = "DualTaskFC")

df_motor_fc <- data.frame(Comp1 = scores_motor_fc[, 1],
                          Comp2 = scores_motor_fc[, 2],
                          Group = group_vector,
                          Block = "SingleMotorFC")

df_behavior <- data.frame(Comp1 = scores_behavior[, 1],
                          Comp2 = scores_behavior[, 2],
                          Group = group_vector,
                          Block = "DualBehavior")

# Combine all into one long dataframe
df_all <- rbind(df_dual_fc, df_motor_fc, df_behavior)

# Set block as a factor with desired order
df_all$Block <- factor(df_all$Block,
                       levels = c("DualTaskFC", "SingleMotorFC", "DualBehavior"))

# Save the data
saveRDS(df_all, "data/df_all_for_Figure7_B.rds")

# Perform Statistical Tests for group comparison
# Run a t-test on Comp1 for each block
stats_results <- df_all %>%
  group_by(Block) %>%
  summarise(
    p_value = t.test(Comp2 ~ Group)$p.value
  ) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
print(stats_results)

# 2. Plot all blocks using facetting with fixed axis
df_all <- readRDS("data/df_all_for_Figure7_B.rds")

p <- ggplot(df_all, aes(x = Comp1, y = Comp2, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.68, type = "norm", linetype = "dashed", size = 0.8) +
  facet_wrap(~ Block, scales = "free") +
  scale_color_manual(values = c("Young" = "#1f77b4", "Old" = "#ff7f0e")) +
  scale_shape_manual(values = c("Young" = 16, "Old" = 17)) +
  labs(
    #title = "Component Scores Across Blocks",
    x = "Component 1",
    y = "Component 2",
    color = "Group",
    shape = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p)

# Save as high-resolution TIFF 
ggsave("20250720_motorNetwork_component_scores_3blocks_final.tiff", 
       plot = p, 
       dpi = 600, 
       width = 10, 
       height = 5, 
       units = "in", 
       compression = "lzw")

# added into supplementary figure

#---------------Splitting the plots per block, full control on each plot
library(patchwork)  # for side-by-side layout

# Split the data into subsets
df_dual <- df_all %>% filter(Block == "DualTaskFC")
df_motor <- df_all %>% filter(Block == "SingleMotorFC")
df_behav <- df_all %>% filter(Block == "DualBehavior")

# 1. DualTaskFC plot
p1 <- ggplot(df_dual, aes(x = Comp1, y = Comp2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.68, linetype = "dashed") +
  coord_cartesian(xlim = c(-12, 25), ylim = c(-15, 15)) +
  labs(title = "Dual Task FC", x = "Component 1", y = "Component 2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 2. SingleMotorFC plot
p2 <- ggplot(df_motor, aes(x = Comp1, y = Comp2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.68, linetype = "dashed") +
  coord_cartesian(xlim = c(-12, 25), ylim = c(-15, 15)) +
  labs(title = "Single Motor FC", x = "Component 1", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 3. DualBehavior plot
p3 <- ggplot(df_behav, aes(x = Comp1, y = Comp2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.68, linetype = "dashed") +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  labs(title = "Dual Behavior", x = "Component 1", y = NULL) +
  theme_minimal(base_size = 14)

# Combine all three plots horizontally with shared legend

final_plot <- (p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "right")

# Print the final combined plot
print(final_plot)

# Save as high-resolution TIFF 
ggsave("20250726_motorNetwork_component_scores_uniformaxis_final.tiff", 
       plot = final_plot, 
       dpi = 600, 
       width = 10, 
       height = 5, 
       units = "in", 
       compression = "lzw")

    
#  Correlation summary panel, version 2, with p value
# Step 1: Data with correlation and p-values
cor_data_all <- data.frame(
  Component = rep(c("Motor Network", "Cognitive Network"), each = 2),
  Block = rep(c("Dual Task FC", "Single Task FC"), times = 2),
  Correlation = c(
    0.432, 0.476,   # Motor
    0.402, 0.465    # Cognitive
  ),
  p_value = c(
    0.0001, 0.0001,    # Motor
    0.0001, 0.0001     # Cognitive
  )
)

# Step 2: Create significance labels based on p-values
cor_data_all <- cor_data_all %>%
  mutate(
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    label = paste0(formatC(Correlation, format = "f", digits = 2), stars)
  )

# Step 3: Factor levels for plot order
cor_data_all$Component <- factor(cor_data_all$Component, levels = c("Motor Network", "Cognitive Network"))
cor_data_all$Block <- factor(cor_data_all$Block, levels = c("Dual Task FC", "Single Task FC"))

# Save the data
saveRDS(cor_data_all,"data/cor_data_all.rds")

# Step 4: Plot heatmap with annotated labels
# Read the data
cor_data_all <- readRDS("data/cor_data_all.rds")

p_corr <- ggplot(cor_data_all, aes(x = Block, y = Component, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 5) +
  scale_fill_gradient2(
    low = "#f7fbff", high = "#08306b", mid = "#c6dbef",
    midpoint = 0.3, limit = c(0.0, 0.6), name = "Correlation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("20250727_brain_behavior_correlation_final.tiff", plot = p_corr, dpi = 600, width = 10, height = 5, units = "in", compression = "lzw")


#---------------------------------------------------
# Plot multiplots for pairwide correlations_ Figure 8
# a simplified version of plotRGCCAMultiPlot():
  plotRGCCAMultiPlot <- function(rgcca_result, group = NULL, block_names = NULL, comp = 1,
                                 lower_fn = "points", upper_fn = "cor", diag_fn = "densityDiag") {
    require(GGally) # Uses standard GGally wrappers ("points", "cor", "densityDiag")
    require(ggplot2)
  
  # Combine variates from each block
  variates_list <- lapply(rgcca_result$variates, function(x) x[, comp])
  df <- as.data.frame(variates_list)
  colnames(df) <- if (is.null(block_names)) names(rgcca_result$variates) else block_names

  # Add group only AFTER selecting columns for ggpairs
  color_mapping <- aes()  # default
  if (!is.null(group)) {
    group <- as.factor(group)
    df$Group <- group
    color_mapping <- aes(color = Group)
    cols_to_plot <- setdiff(names(df), "Group")
  } else {
    cols_to_plot <- names(df)
  }
  
  # ggpairs plot
  p <-GGally::ggpairs(
    df,
    columns = cols_to_plot,
    mapping = color_mapping,
    lower = list(continuous = wrap(lower_fn, alpha = 0.7, size = 2)),
    upper = list(continuous = wrap(upper_fn, size = 4)),
    diag  = list(continuous = wrap(diag_fn, alpha = 0.5))
  )
  
  print(p)
  return(p)
  }
  
  # call the plotting function for pairwise correlations_ Figure 8
  myplot <- plotRGCCAMultiPlot(
    rgcca_result = res_rgcca_dual,
    group = group_vector,  
    block_names = c("DualTaskFC", "SingleMotorFC", "DualBehavior"),
    comp = 1
    #diag_fn = diag_no_grid
  )

  # Save as high-resolution TIFF (good for journals)
  ggsave("20250722_motorNetwork_correlationMatrix_comp2_final.tiff", plot = myplot, dpi = 600, width = 10, height = 5, units = "in", compression = "lzw")
  
  # Plot both comp 1 and 2.  And comp 1 is much better to distinguish age group.
  
  
# - plot loadings of behavior. ------------------------------------
# Reusable Function from rgcca loadings
install.packages("tibble") # If haven’t installed it yet
  
# Function to plot loadings from RGCCA results: Figure 9_Behavior loading
plot_behavior_loadings <- function(rgcca_result, comp = 1, block_name = "DualBehavior", motor_keywords = c("Motor", "motor")) {
  
  # Extract loadings for the specified block and component
  loadings <- rgcca_result$loadings[[block_name]][, comp, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Variable")
  
  colnames(loadings)[2] <- "Loading"
  
  # Add category based on variable name
  loadings$Category <- ifelse(
    grepl(paste(motor_keywords, collapse = "|"), loadings$Variable),
    "Motor",
    "Cognitive"
  )
  
  # Sort variables by absolute loading
  loadings <- loadings %>%
    arrange(desc(abs(Loading))) %>%
    mutate(Variable = factor(Variable, levels = rev(Variable)))  # For proper ggplot ordering
  
  # Plot
  p <- ggplot(loadings, aes(x = Variable, y = Loading, fill = Category)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Motor" = "#1f77b4", "Cognitive" = "#2ca02c")) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Behavioral Loadings (Component", comp, ")"),
      x = NULL,
      y = "Loading Weight",
      fill = "Task Type"
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  print(p)
}

# Plot behavior loading
my_plot <- plot_behavior_loadings(res_rgcca_dual, comp = 1)

# Save as high-res TIFF
ggsave("2025_0713_Behavior_Loadings_Comp1_sort.tiff", 
       plot = my_plot, 
       width = 8, 
       height = 5, 
       dpi = 600, 
       device = "tiff", 
       compression = "lzw")

# -------------------------------------------------------
# Pyramid Barplot Function for the daultask FC

plot_pyramid_loadings <- function(rgcca_result, block_name, comp = 2, top_n = 20) {
  # Extract and clean loadings
  loadings <- rgcca_result$loadings[[block_name]][, comp, drop = FALSE] %>%
    as.data.frame()
  loadings$FullName <- rownames(loadings)
  colnames(loadings)[1] <- "Loading"
  
  pattern_suffix <- if (block_name == "SingleMotorFC") "Motor_single" else "Go_dual"
  # Clean variable labels
  loadings$Variable <- loadings$FullName %>%
    str_extract(paste0("M_ROI_.*?\\.and\\.M_ROI_.*?\\.at\\.", pattern_suffix)) %>%
    str_replace_all("M_ROI_", "") %>%
    str_replace_all("\\.and\\.", " – ") %>%
    str_replace(paste0("\\.at\\.", pattern_suffix), "")
  
  # Keep top N by absolute value
  top_loadings <- loadings %>%
    slice_max(order_by = abs(Loading), n = top_n) %>%
    arrange(Loading) %>%
    mutate(Variable = factor(Variable, levels = Variable))  # Keep proper Y ordering
  
  # Plot mirrored barplot
  p <- ggplot(top_loadings, aes(x = Loading, y = Variable, fill = Loading > 0)) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = c("TRUE" = "#66c2a5", "FALSE" = "#8da0cb"),  # softer teal and blue
      guide = "none") +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    labs(
      title = paste("Barplot of Loadings -", block_name, "(Component", comp, ")"),
      x = "Loading Weight",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  print(p)
  return(p)
}

my_plot <- plot_pyramid_loadings(res_rgcca_dual, block_name = "DualTaskFC", comp = 1)
my_plot <- plot_pyramid_loadings(res_rgcca_dual, block_name = "SingleMotorFC", comp = 1)

ggsave("20250719_DualTaskFC_loadings_motor_network_comp1_600.tiff", 
       my_plot, width = 10, height = 6, dpi = 600, compression = "lzw")

# __ End __#