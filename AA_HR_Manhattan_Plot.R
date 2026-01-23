# Manhattan plots for Hazard Ratios by Amino Acid Position and Locus
# Adapted script for CSV data with columns: Locus, Position, Don_rec_AA, HR, Conf_int_95, Var_p, Model_p, N, Diag_meas

# Load required packages
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}
if (!require("ggrepel")) {
  install.packages("ggrepel")
  library("ggrepel")
}
if (!require("ggtext")) {
  install.packages("ggtext")
  library("ggtext")
}

data <- read.csv("HR by don rec AAMM.csv", stringsAsFactors = FALSE)
check_data_types <- function(df) {
  cat("Data types in the input file:\n")
  cat("Locus:", class(df$Locus), "\n")
  cat("Position:", class(df$Position), "\n")
  cat("Don_rec_AA:", class(df$Don_rec_AA), "\n")
  cat("HR:", class(df$HR), "\n")
  cat("Conf_int_95:", class(df$Conf_int_95), "\n")
  cat("Var_p:", class(df$Var_p), "\n")
  cat("Model_p:", class(df$Model_p), "\n")
  cat("N:", class(df$N), "\n")
  cat("Diag_meas:", class(df$Diag_meas), "\n")
}
check_data_types(data)

# Function to handle special p-value formats like "<.0001", "<0.001", etc.
parse_pvalue <- function(pvalue_string) {
  pvalue_string <- as.character(pvalue_string)
  
  if (is.na(pvalue_string) || pvalue_string == "" || pvalue_string == "NA") {
    return(NA)
  }
  
  # Handle formats like "<.0001", "<0.001", "< 0.001", etc.
  if (grepl("^\\s*<\\s*\\.?0*1*e?-?[0-9]*$", pvalue_string, ignore.case = TRUE)) {
    # Extract the number after < symbol
    number_part <- gsub("^\\s*<\\s*", "", pvalue_string)
    
    # Handle cases like ".0001" -> "0.0001"
    if (startsWith(number_part, ".")) {
      number_part <- paste0("0", number_part)
    }
    
    # Convert to numeric and make it slightly smaller than the threshold
    threshold <- as.numeric(number_part)
    if (!is.na(threshold)) {
      return(threshold * 0.99)
    }
  }
  
  numeric_val <- suppressWarnings(as.numeric(pvalue_string))
  return(numeric_val)
}

data$Position <- as.numeric(data$Position)
data$HR <- as.numeric(data$HR)
data$Var_p <- sapply(data$Var_p, parse_pvalue)
data$Model_p <- sapply(data$Model_p, parse_pvalue)
data$N <- as.integer(data$N)
data$Diag_meas <- as.numeric(data$Diag_meas)

# Show summary of converted p-values
cat("Summary of Var_p after conversion:\n")
cat("Range:", range(data$Var_p, na.rm = TRUE), "\n")
cat("Number of values < 0.05:", sum(data$Var_p < 0.05, na.rm = TRUE), "\n")
cat("Number of NAs:", sum(is.na(data$Var_p)), "\n")

cat("\nSummary of Model_p after conversion:\n")
cat("Range:", range(data$Model_p, na.rm = TRUE), "\n")
cat("Number of values < 0.05:", sum(data$Model_p < 0.05, na.rm = TRUE), "\n")
cat("Number of NAs:", sum(is.na(data$Model_p)), "\n")

check_data_types(data)

# Function to parse confidence intervals and extract lower and upper bounds
parse_confidence_interval <- function(ci_string) {
  # Extract numbers from string like "(0.97-1.04)"
  ci_clean <- gsub("[()]", "", ci_string)  # Remove parentheses
  ci_parts <- strsplit(ci_clean, "-")[[1]]
  return(list(
    lower = as.numeric(ci_parts[1]),
    upper = as.numeric(ci_parts[2])
  ))
}

data$CI_lower <- sapply(data$Conf_int_95, function(x) parse_confidence_interval(x)$lower)
data$CI_upper <- sapply(data$Conf_int_95, function(x) parse_confidence_interval(x)$upper)

manhattan_plot_hr <- function(locus_data, locus_name) {
  
  unique_positions <- sort(unique(locus_data$Position))
  position_mapping <- data.frame(
    Position = unique_positions,
    x_position = 1:length(unique_positions)
  )
  
  plot_data <- locus_data %>%
    left_join(position_mapping, by = "Position") %>%
    arrange(Position) %>%
    mutate(
      base_label = paste0(Don_rec_AA, " (N=", N, ")"),
      var_p_sig = ifelse(Var_p < 0.05, "*", ""),
      label_with_p = paste0(base_label, var_p_sig),
      is_highlighted = Diag_meas > 2,
      model_p_sig = ifelse(Model_p < 0.05, "*", ""),
      position_label = paste0(Position, model_p_sig),
      ci_excludes_1 = (CI_lower > 1) | (CI_upper < 1),
      final_label = ifelse(ci_excludes_1, label_with_p, "")
    )
  
  position_labels <- plot_data %>%
    group_by(Position, x_position) %>%
    summarise(
      model_p_sig = first(model_p_sig),
      position_label = first(position_label),
      .groups = "drop"
    ) %>%
    arrange(x_position)
  
  y_min <- min(c(plot_data$HR, plot_data$CI_lower), na.rm = TRUE)
  y_max <- max(c(plot_data$HR, plot_data$CI_upper), na.rm = TRUE)
  y_range_padding <- (y_max - y_min) * 0.1
  y_min_padded <- max(0.1, y_min - y_range_padding)  # Don't go below 0.1 for log scale
  y_max_padded <- y_max + y_range_padding
  
  p <- ggplot(plot_data, aes(x = x_position, y = HR)) +
    
    geom_point(aes(color = is_highlighted), 
               size = 2, alpha = 0.8) +
    
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
    
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
    
    geom_text_repel(aes(label = final_label,
                        color = is_highlighted),
                    size = 3,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    max.overlaps = Inf,
                    show.legend = FALSE) +
    
    scale_x_continuous(
      breaks = position_labels$Position,
      labels = NULL
    ) +
    
    scale_y_continuous(
      name = "Hazard Ratio (95% CI)",
      limits = c(y_min_padded, y_max_padded)
    ) +
    
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    
    ggtitle(paste("Hazard Ratios by Amino Acid Position - Locus", locus_name)) +
    xlab("Position") +

  return(p)
}

unique_loci <- unique(data$Locus)
for (locus in unique_loci) {
  locus_data <- data %>% filter(Locus == locus)
  
  cat("Creating plot for locus:", locus, "with", nrow(locus_data), "data points\n")
  
  plot <- manhattan_plot_hr(locus_data, locus)
  filename <- paste0("Don_Rec_AA_MM_permutation_HR_manhattan_", locus, ".png")
  ggsave(filename, plot, width = 12, height = 8, dpi = 300)
  
  cat("Saved plot as:", filename, "\n")
}

create_combined_plot <- function(data) {
  data$lgroup <- 'NA'
  data[data$Locus == 'A', "lgroup"] <- 'L1'
  data[data$Locus == 'C', "lgroup"] <- 'L2'
  data[data$Locus == 'B', "lgroup"] <- 'L3'
  data[data$Locus == 'DRB345', "lgroup"] <- 'L4'
  data[data$Locus == 'DRB1', "lgroup"] <- 'L5'
  data[data$Locus == 'DQA1', "lgroup"] <- 'L6'
  data[data$Locus == 'DQB1', "lgroup"] <- 'L7'
  
  combined_data <- data %>%
    group_by(lgroup) %>%
    summarise(max_pos = max(Position, na.rm = TRUE), .groups = "drop") %>%
    mutate(cumsum_pos = cumsum(c(0, max_pos[-length(max_pos)]))) %>%
    right_join(data, by = "lgroup") %>%
    mutate(
      cumulative_position = Position + cumsum_pos,
      base_label = paste0(Don_rec_AA, " (N=", N, ")"),
      var_p_sig = ifelse(Var_p < 0.05, "*", ""),
      label_with_p = paste0(base_label, var_p_sig),
      is_highlighted = Diag_meas > 2
    )
  
  axis_labels <- combined_data %>%
    group_by(Locus) %>%
    summarise(center = mean(cumulative_position), .groups = "drop")
  
  p <- ggplot(combined_data, aes(x = cumulative_position, y = HR)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  width = 0.5, alpha = 0.6, color = "gray60") +
    geom_point(aes(color = as.factor(lgroup)), size = 2, alpha = 0.8) +
    scale_color_manual(values = rep(c("blue", "red"), 4)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_x_continuous(
      breaks = axis_labels$center,
      labels = axis_labels$Locus
    ) +
    scale_y_continuous(
      name = "Hazard Ratio (95% CI)",
      breaks = c(0.5, 0.7, 1.0, 1.5, 2.0, 3.0)
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    ggtitle("Combined Manhattan Plot - All Loci") +
    xlab("Locus")
  
  return(p)
}

combined_plot <- create_combined_plot(data)
ggsave("Don_Rec_AA_MM_permutation_HR_manhattan_ALL.png", combined_plot, width = 16, height = 8, dpi = 300)
