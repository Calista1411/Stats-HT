# This function plots the single nucleotide counts across windows of a specified
# length, for sequence files within two folders. Allows us to compare the base 
# composition of descendant sequences generated from different ancestral sequences
# Arguments:
  # biased_dir: path to folder containing descendant sequence files for the non-random
              # (e.g. coding DNA) ancestral sequences
  # random_dir: path to folder containing descendant sequence files for the randomly
              # sampled ancestral sequences
  # window: window size (in bases) used to aggregate counts, defaults to 1000
  # pattern: file format, defaults to text files
# Output: Line plots showing the counts of base per window for each sequence
        # file, along with the corresponding mean line for the sequences in each folder
        # Separate line plots are produced for each base (A, C, G, T)

plot_single_nuc_counts <- function(
    biased_dir,
    random_dir,
    window  = 1000,
    pattern = "\\.txt$"
) {

  files_b <- list.files(biased_dir, pattern = pattern, full.names = TRUE)
  files_r <- list.files(random_dir, pattern = pattern, full.names = TRUE)
  
  # Helper  function to split sequence into windows and count the numbers of A,C,G,T
  make_windows <- function(n, window) split(seq_len(n), ceiling(seq_len(n) / window))
  count_mono <- function(chars_vec) {
    tab <- table(factor(chars_vec, levels = 1:4))
    as.integer(tab)
  }
  
  summarise_windows <- function(num_vec, file_id, group, window) {
    wnd <- make_windows(length(num_vec), window)
    purrr::map_dfr(seq_along(wnd), function(wi) {
      seg <- num_vec[wnd[[wi]]]
      cnt <- count_mono(seg)
      tibble(
        file = file_id,
        group = group,
        window = wi,
        base = c("A","C","G","T"),
        count = cnt
      )
    })
  }
  
  read_group <- function(paths, group_label) {
    purrr::imap_dfr(paths, function(pth, i) {
      seq_num <- read_seq(pth)
      max_len <- floor(length(seq_num) / window) * window
      if (max_len > 0 && max_len < length(seq_num)) {
        seq_num <- seq_num[seq_len(max_len)]
      }
      fid <- sprintf("%s-%02d", group_label, i)
      summarise_windows(seq_num, fid, group_label, window)
    })
  }
  
  df <- bind_rows(
    read_group(files_b, "Biased"),
    read_group(files_r, "Random")
  )
  
  # Compute mean counts per base, window, and group
  df_mean <- df |>
    group_by(group, base, window) |>
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

  colors <- c(
    "Genome replicates" = "#f4a6a6",
    "Random replicates" = "#9ecae1",
    "Genome mean"       = "#b22222",
    "Random mean"       = "#1f77b4"
  )
  
  make_plot_for_base <- function(b) {
    d_rep <- df %>% filter(base == b)
    d_m   <- df_mean %>% filter(base == b)
    
    ggplot() +
      geom_line(data = d_rep %>% filter(group == "Biased"),
                aes(x = window, y = count, group = file, color = "Genome replicates"),
                alpha = 0.5, linewidth = 0.5) +
      geom_line(data = d_rep %>% filter(group == "Random"),
                aes(x = window, y = count, group = file, color = "Random replicates"),
                alpha = 0.5, linewidth = 0.5) +
      geom_line(data = d_m %>% filter(group == "Biased"),
                aes(x = window, y = mean_count, color = "Genome mean"),
                linewidth = 1.2) +
      geom_line(data = d_m %>% filter(group == "Random"),
                aes(x = window, y = mean_count, color = "Random mean"),
                linewidth = 1.2) +
      scale_color_manual(values = colors) +
      labs(
        title = sprintf("%s counts per 1000 bases", b),
        x = "Window index",
        y = "Count per 1000 bases",
        color = NULL
      ) +
      theme_bw(base_family = "CMU Serif", base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
  }
  
  list(
    A_plot = make_plot_for_base("A"),
    C_plot = make_plot_for_base("C"),
    G_plot = make_plot_for_base("G"),
    T_plot = make_plot_for_base("T")
  )
}