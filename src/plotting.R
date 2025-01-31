library(dplyr)
library(tidyr)
library(ggplot2)


plot_metrics_by_param <- function(df, param_col, metric_type, log_scale = FALSE, legend_title = "Metric") {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Filter columns for the chosen metric type
  metric_cols <- names(df)[grepl(metric_type, names(df), ignore.case = TRUE)]
  
  # Remove columns with all NA values
  metric_cols <- metric_cols[sapply(df[metric_cols], function(col) !all(is.na(col)))]
  
  # Group by the parameter column and calculate mean and standard deviation
  grouped_df <- df %>%
    group_by_at(param_col) %>%
    summarise(across(
      all_of(metric_cols), 
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    )) %>%
    ungroup()
  
  # Reshape the data for plotting
  mean_data <- grouped_df %>%
    pivot_longer(cols = ends_with("_mean"), 
                 names_to = "Metric", 
                 values_to = "Value",
                 names_pattern = "(.*)_mean")
  
  std_data <- grouped_df %>%
    pivot_longer(cols = ends_with("_sd"), 
                 names_to = "Metric", 
                 values_to = "std",
                 names_pattern = "(.*)_sd")
  
  # Merge mean and standard deviation data
  plot_data <- left_join(mean_data, std_data, by = c(param_col, "Metric"))
  
  # print(plot_data)
  
  # Set custom y-axis label based on metric type
  y_label <- ifelse(grepl("sensitivity", metric_type, ignore.case = TRUE),
                    "Recall",
                    ifelse(grepl("precision", metric_type, ignore.case = TRUE),
                           "Precision",
                           paste(metric_type, "Value")))
  
  # Define colors for both Recall and Precision metrics
  color_map <- c(
    "Precision_Baseline_25" = "green",
    "Precision_CV" = "red",
    "Precision_Baseline_50" = "blue",
    "Precision_Baseline_75" = "orange",
    "Sensitivity_CV" = "red",
    "Sensitivity_Baseline_25" = "green",
    "Sensitivity_Baseline_50" = "blue",
    "Sensitivity_Baseline_75" = "orange"
  )
  

  p <- ggplot(plot_data, aes_string(x = param_col, y = "Value", color = "Metric", group = "Metric")) +
    geom_line(size=3) +
    geom_point(size=4) +
    geom_errorbar(aes(ymin = pmax(0, Value - std), ymax = pmin(1, Value + std)), width = 0.08, size=1) + # Clipped error bars
    #   scale_y_continuous(limits = c(0, 1)) +  # Enforce same y-axis scale
    scale_color_manual(values = color_map) + 
    labs(
      x = param_col,
      y = y_label,
      color = legend_title
    ) +
    theme_minimal() +
    ylim(0,1) + 
    theme(
      axis.title = element_text(face = "bold", size = 25),
      axis.text = element_text(face = "bold", size = 18),
 
      panel.grid.major = element_line(color = "gray70", size = 0.5),  # Darker and thicker major grid
      panel.grid.minor = element_line(color = "gray70", size = 0.5),  
   #   panel.border = element_rect(color = "gray10", fill = NA, size = 1)
    )
  
  # Apply log scale if specified
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  
  # Print the plot
  print(p)
}



plot_metrics_by_param_snr <- function(df, param_col, metric_type, log_scale = FALSE, legend_title = "Metric", custom_x_labels = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Filter columns for the chosen metric type
  metric_cols <- names(df)[grepl(metric_type, names(df), ignore.case = TRUE)]
  
  # Remove columns with all NA values
  metric_cols <- metric_cols[sapply(df[metric_cols], function(col) !all(is.na(col)))]
  
  # Group by the parameter column and calculate mean and standard deviation
  grouped_df <- df %>%
    group_by_at(param_col) %>%
    summarise(across(
      all_of(metric_cols), 
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    )) %>%
    ungroup()
  
  # Reshape the data for plotting
  mean_data <- grouped_df %>%
    pivot_longer(cols = ends_with("_mean"), 
                 names_to = "Metric", 
                 values_to = "Value",
                 names_pattern = "(.*)_mean")
  
  std_data <- grouped_df %>%
    pivot_longer(cols = ends_with("_sd"), 
                 names_to = "Metric", 
                 values_to = "std",
                 names_pattern = "(.*)_sd")
  
  # Merge mean and standard deviation data
  plot_data <- left_join(mean_data, std_data, by = c(param_col, "Metric"))
  
  # Set custom y-axis label based on metric type
  y_label <- ifelse(grepl("sensitivity", metric_type, ignore.case = TRUE),
                    "Recall",
                    ifelse(grepl("precision", metric_type, ignore.case = TRUE),
                           "Precision",
                           paste(metric_type, "Value")))
  
  # Define colors for both Recall and Precision metrics
  color_map <- c(
    "Precision_Baseline_25" = "green",
    "Precision_CV" = "red",
    "Precision_Baseline_50" = "blue",
    "Precision_Baseline_75" = "orange",
    "Sensitivity_CV" = "red",
    "Sensitivity_Baseline_25" = "green",
    "Sensitivity_Baseline_50" = "blue",
    "Sensitivity_Baseline_75" = "orange"
  )
  
  # Convert x-axis column to character to ensure discrete axis
  plot_data[[param_col]] <- as.character(plot_data[[param_col]])
  
  # Create the base plot
  p <- ggplot(plot_data, aes_string(x = param_col, y = "Value", color = "Metric", group = "Metric")) +
    geom_line(size=3) +
    geom_point(size=4) +
    geom_errorbar(aes(ymin = pmax(0, Value - std), ymax = pmin(1, Value + std)), width = 0.25, size=1) + # Clipped error bars
    scale_color_manual(values = color_map) + 
    labs(
      x = "SNR",
      y = y_label,
      color = legend_title
    ) +
    theme_minimal() +
    ylim(0, 1) + 
    theme(
      axis.title = element_text(face = "bold", size = 25),
      axis.text = element_text(face = "bold", size = 18),
      #  legend.text = element_text(face = "bold", size = 16),
      # legend.title = element_text( face = "bold", size = 18),
      #  plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      # panel.grid.major = element_line(size = 1),
      # panel.grid.minor = element_blank()  # Remove minor grid lines for clarity
      panel.grid.major = element_line(color = "gray70", size = 0.5),  # Darker and thicker major grid
      panel.grid.minor = element_line(color = "gray70", size = 0.5),  
      #   panel.border = element_rect(color = "gray10", fill = NA, size = 1)
      legend.position = "bottom",      # Move legend below the plot
      legend.direction = "horizontal"
    )
  
  # Apply custom x-axis labels if provided
  if (!is.null(custom_x_labels)) {
    p <- p + scale_x_discrete(labels = custom_x_labels)
  }
  
  # Apply log scale if specified (Note: Log scale doesn't work with discrete axes)
  if (log_scale) {
    message("Log scale is disabled because x-axis is treated as discrete.")
  }
  
  # Print the plot
  print(p)
}


plot_metrics_by_param_rho_sweep <- function(df, param_col, metric_type, log_scale = FALSE, legend_title = "Metric") {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Filter columns for the chosen metric type
  metric_cols <- names(df)[grepl(metric_type, names(df), ignore.case = TRUE)]
  
  # Remove columns with all NA values
  metric_cols <- metric_cols[sapply(df[metric_cols], function(col) !all(is.na(col)))]
  
  # Group by the parameter column and calculate mean and standard deviation
  grouped_df <- df %>%
    group_by_at(param_col) %>%
    summarise(across(
      all_of(metric_cols), 
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    )) %>%
    ungroup()
  
  # Reshape the data for plotting
  mean_data <- grouped_df %>%
    pivot_longer(cols = ends_with("_mean"), 
                 names_to = "Metric", 
                 values_to = "Value",
                 names_pattern = "(.*)_mean")
  
  std_data <- grouped_df %>%
    pivot_longer(cols = ends_with("_sd"), 
                 names_to = "Metric", 
                 values_to = "std",
                 names_pattern = "(.*)_sd")
  
  # Merge mean and standard deviation data
  plot_data <- left_join(mean_data, std_data, by = c(param_col, "Metric"))
  
  # print(plot_data)
  
  # Set custom y-axis label based on metric type
  y_label <- ifelse(grepl("Prop", metric_type, ignore.case = TRUE),
                    "Proportion ϕ",
                    ifelse(grepl("Perc", metric_type, ignore.case = TRUE),
                           "Percent better",
                           paste(metric_type, "Value")))
  
  # Define colors for both Recall and Precision metrics
  color_map <- c(
    "Prop_ss" = "red",
    "Prop_25" = "green",
    "Prop_50" = "blue",
    "Prop_75" = "orange",
    "Perc_ss" = "red",
    "Perc_25" = "green",
    "Perc_50" = "blue",
    "Perc_75" = "orange"
  )
  

  if (metric_type=="Prop"){
    p <- ggplot(plot_data, aes_string(x = param_col, y = "Value", color = "Metric", group = "Metric")) +
      geom_line(size=3) +
      geom_point(size=4) +
      geom_errorbar(aes(ymin = pmax(0, Value - std), ymax = pmin(1, Value + std)), width = 0.02, size=1) + # Clipped error bars
      #   scale_y_continuous(limits = c(0, 1)) +  # Enforce same y-axis scale
      scale_color_manual(values = color_map) + 
      labs(
        x = param_col,
        y = y_label,
        color = legend_title
      ) +
      theme_minimal() +
      ylim(0,1) + 
      theme(
        axis.title = element_text(face = "bold", size = 25),
        axis.text = element_text(face = "bold", size = 18),
        #  legend.text = element_text(face = "bold", size = 16),
        # legend.title = element_text( face = "bold", size = 18),
        #  plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
        # panel.grid.major = element_line(size = 1),
        # panel.grid.minor = element_blank()  # Remove minor grid lines for clarity
        panel.grid.major = element_line(color = "gray70", size = 0.5),  # Darker and thicker major grid
        panel.grid.minor = element_line(color = "gray70", size = 0.5),  
        #   panel.border = element_rect(color = "gray10", fill = NA, size = 1)
        # legend.position = "bottom",      # Move legend below the plot
        # legend.direction = "horizontal"
      )
    
  } else {
  }

  
  # Apply log scale if specified
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  
  # Print the plot
  print(p)
}



