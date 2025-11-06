# This function plots the AIC and BIC values for Markov chain models of
# different orders, across sequences of varying lengths
# Arguments: 
  # ... : one or more data frames containing AIC and BIC results for a collection 
        # of sequences
  # labels: sequence length labels corresponding to the input dataframes
  # order: integer vector specifying which Markov orders to plot
  # colours: colour codes for AIC and BIC lines 
  # base_family: base font family used for text elements, defaults to "CMU Serif"
  # line_alpha: transparency level for plotted lines, defaults to 0.45
  # line_width: line width for AIC and BIC lines, defaults to 0.7
  # point_size: size of scatter points plotted along each line, defaults to 0.8
# Output: AIC and BIC plots for all provided datasets

plot_ic <- function(
    ...,
    labels,
    orders      = 0:5,
    colours     = c(AIC = "#FFB6C1", BIC = "#D3B5E5"),
    base_family = "CMU Serif",
    line_alpha  = 0.45,
    line_width  = 0.7,
    point_size  = 0.8
) {
  
  results_list <- list(...)
  res_all <- bind_rows(
    map2(results_list, labels, \(x, lab) mutate(x, length = lab))
  )
  
  # Convert to long format
  long_ic <- res_all |>
    select(any_of(c("file", "length")), matches("^(AIC|BIC)\\d+$")) |>
    pivot_longer(
      cols = matches("^(AIC|BIC)\\d+$"),
      names_to = c("criterion", "order"),
      names_pattern = "(AIC|BIC)(\\d+)",
      values_to = "value"
    ) |>
    mutate(
      order   = as.integer(order),
      value   = suppressWarnings(as.numeric(value)),
      file_id = paste(length, file, sep = "::")
    ) |>
    filter(order %in% orders) |>
    arrange(file_id, criterion, order)
  
  # Build one plot panel per sequence length
  build_panel <- function(len_label) {
    long_ic |>
      filter(length == len_label) |>
      ggplot(aes(x = order, y = value,
                 group = interaction(file_id, criterion),
                 colour = criterion)) +
      geom_line(alpha = line_alpha, linewidth = line_width) +
      geom_point(alpha = line_alpha, size = point_size) +
      facet_wrap(~ criterion, nrow = 1, scales = "fixed") +
      scale_x_continuous(breaks = sort(unique(orders))) +
      scale_colour_manual(values = colours, guide = "none") +
      labs(title = len_label, x = "Markov Order", y = "AIC / BIC Value") +
      theme_minimal(base_size = 12, base_family = base_family) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, face = "bold",
                                  family = base_family, size = 13,
                                  margin = margin(b = 6)),
        axis.title.x = element_text(margin = margin(t = 14), family = base_family),
        axis.title.y = element_text(margin = margin(r = 14), family = base_family),
        strip.text = element_text(face = "bold", family = base_family)
      )
  }
  
  panels <- map(labels, build_panel)
  wrap_plots(panels, ncol = 1, heights = rep(1, length(labels)))
}