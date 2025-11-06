# This function plots percentage of bases that differ between the ancestral and 
# descendant sequence across different evolution times
# Arguments: 
  # file: path to CSV file containing a column identifying replicates and one 
        # column per evolution time, with each cell having the % of bases changed
# Output: 
  # Line graph of % of bases changed against evolution time, with light pink lines
  # showing the trend for individual replicates, and purple line and points showing
  # mean % across replicates

plot_bases_changed <- function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  
  df_long <- df |>
    pivot_longer(
      cols = matches("^t="),
      names_to  = "time",
      values_to = "percent_changed"
    ) |>
    mutate(
      time = readr::parse_number(time),
      percent_changed = as.numeric(percent_changed)
    )
  
  # Calculate mean across replicates
  df_mean <- df_long |>
    group_by(time) |>
    summarise(mean_change = mean(percent_changed, na.rm = TRUE), .groups = "drop")
  
  # Plot line graph (light pink lines for each replicate, purple line for mean)
  ggplot(df_long, aes(x = time, y = percent_changed, group = replicate_id)) +
    geom_line(color = "#F8BBD0", linewidth = 0.7, alpha = 0.7) +   
    geom_line(data = df_mean, aes(x = time, y = mean_change),
              inherit.aes = FALSE, color = "#9C27B0", linewidth = 1.4) +  
    geom_point(data = df_mean, aes(x = time, y = mean_change),
               inherit.aes = FALSE, color = "#9C27B0", size = 2.5) +
    scale_x_continuous(breaks = sort(unique(df_long$time))) +
    labs(
      x = "Evolution time (t)",
      y = "% of bases changed"
    ) +
    theme_minimal(base_family = "CMU Serif", base_size = 14) +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text  = element_text(size = 11),
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}