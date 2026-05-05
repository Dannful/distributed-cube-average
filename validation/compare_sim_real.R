library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(optparse)

option_list <- list(
  make_option("--real-data-dir", type = "character", default = NULL),
  make_option("--simulation-data-dir", type = "character", default = NULL),
  make_option("--problem-sizes", type = "character", default = NULL),
  make_option("--time-max", type = "integer", default = NULL)
)
opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$`real-data-dir`) && is.null(opts$`simulation-data-dir`)) {
  stop("At least one of --real-data-dir or --simulation-data-dir must be provided.")
}
if (is.null(opts$`problem-sizes`)) {
  stop("--problem-sizes is required.")
}

real_traces_path <- opts$`real-data-dir`
sim_traces_path <- opts$`simulation-data-dir`
problem_sizes <- as.numeric(str_split_1(opts$`problem-sizes`, ","))
has_real <- !is.null(real_traces_path)
has_sim <- !is.null(sim_traces_path)
time_max <- opts$`time-max`

read_dataset <- function(file_path) {
  dta <- read_csv(file_path, show_col_types = FALSE, col_names = c(
    "Nature", "Container",
    "Type", "Start", "End", "Duration", "Imbrication", "Value"
  ))
  dta |>
    select(-Type, -Imbrication, -Nature) |>
    mutate(Container = as.integer(gsub("rank-?", "", Container)), Value = gsub(
      "^PMPI_",
      "MPI_", Value
    )) |>
    rename(Rank = Container, Operation = Value) -> df.states

  df.states
}

better_contrast_style <- function() {
  scale_fill_manual(
    values = c(
      MPI_Irecv = "#4DAF4A",
      MPI_Isend = "#E41A1C",
      MPI_Waitall = "#377EB8",
      Compute = "#FFD92F"
    )
  )
}

load_simulation <- function(traces_path) {
  csv_files <- list.files(
    recursive = TRUE, full.names = TRUE, pattern = "\\.csv$",
    path = traces_path
  )
  df <- csv_files |>
    set_names() |>
    map_dfr(read_dataset, .id = "source_file") |>
    mutate(problem_size = as.integer(str_extract(source_file, "\\d+"))) |>
    select(-source_file)
  df
}

load_real_traces <- function(traces_path) {
  rst_files <- list.files(
    recursive = TRUE, full.names = TRUE, pattern = "\\.rst$",
    path = traces_path
  )
  rst_directories <- dirname(rst_files)
  rst_grouped <- split(rst_files, rst_directories)
  df <- lapply(names(rst_grouped), \(rst_directory) {
    name_regex <- "\\d+-\\d+-\\d+/{1,2}(\\d+)/{1,2}(\\d+)"
    matches <- str_match(rst_directory, name_regex)
    problem_size <- as.integer(matches[, 3])
    run_number <- as.integer(matches[, 2])
    all_rst_files <- paste0(rst_directory, "/*.rst")
    system(paste("aky_converter", "-l", all_rst_files, ">", "dc.trace"))
    system(paste(
      "pj_dump", "-z", "-l", "9", "dc.trace", "|", "grep", "^State",
      ">", "dc.csv"
    ))
    read_dataset("dc.csv") |>
      mutate(run_number = run_number) |>
      mutate(problem_size = problem_size)
  })
  file.remove("dc.trace")
  file.remove("dc.csv")
  bind_rows(df)
}

generate_gantt_chart <- function(df) {
  df |>
    ggplot() +
    geom_rect(aes(
      xmin = Start, xmax = End, ymin = Rank - 0.5,
      ymax = Rank + 0.5, fill = Operation
    )) +
    scale_y_continuous(breaks = scales::breaks_width(1)) +
    xlab("Time [seconds]") +
    ylab("Rank") +
    theme_bw(base_size = 14) +
    theme(
      plot.margin = unit(c(10, 25, 10, 10), "pt"),
      legend.margin = margin(t = 0, unit = "cm"), panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left", legend.box.spacing = unit(
        0,
        "pt"
      ), legend.box.margin = margin(0, 0, 0, 0), legend.title = element_text(size = 10)
    ) +
    better_contrast_style()
}

generate_combined_gantt_chart <- function(real_df, sim_df) {
  bind_rows(list(real = real_df, sim = sim_df), .id = "src") |>
    mutate(src = as.factor(src)) |>
    generate_gantt_chart() +
    facet_wrap(~src, labeller = as_labeller(c("real" = "Real", "sim" = "Simulation")))
}

compute_load_data <- function(df, n_slices = 100) {
  t_min <- min(df$Start)
  t_max <- max(df$End)
  slice_width <- (t_max - t_min) / n_slices
  slices <- tibble(
    Slice = seq_len(n_slices),
    StartTS = t_min + (seq_len(n_slices) - 1) * slice_width,
    EndTS = t_min + seq_len(n_slices) * slice_width,
    DurationTS = slice_width
  )
  df |>
    cross_join(slices) |>
    mutate(Time.Sum = pmax(0, pmin(End, EndTS) - pmax(Start, StartTS))) |>
    filter(Time.Sum > 0) |>
    group_by(Rank, Slice, StartTS, EndTS, DurationTS, Operation) |>
    summarise(Time.Sum = sum(Time.Sum), .groups = "drop")
}

generate_load_chart <- function(df) {
  load_theme <- list(
    theme_bw(base_size = 14),
    theme(
      panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_text(size = 10)
    ),
    better_contrast_style()
  )
  yconfm <- df |>
    distinct(Rank) |>
    arrange(Rank) |>
    mutate(Position = 0:(n() - 1))
  p_rank <- df |>
    filter(Time.Sum != 0) |>
    group_by(Rank, Slice, StartTS, EndTS, DurationTS, Operation) |>
    summarize(Time.Sum = sum(Time.Sum), .groups = "drop_last") |>
    ungroup(Operation) |>
    arrange(Operation) |>
    mutate(TaskHeight = Time.Sum / DurationTS) |>
    mutate(TaskPosition = cumsum(TaskHeight) - TaskHeight) |>
    left_join(yconfm, by = "Rank") |>
    ggplot() +
    scale_y_continuous(breaks = yconfm$Position, labels = yconfm$Rank) +
    geom_rect(alpha = 0.5, aes(
      fill = Operation,
      xmin = StartTS, xmax = EndTS,
      ymin = Position + TaskPosition,
      ymax = Position + TaskPosition + TaskHeight
    )) +
    ylab("Rank [index]") +
    load_theme +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p_load <- df |>
    filter(Time.Sum != 0) |>
    group_by(Slice, StartTS, EndTS, DurationTS, Operation) |>
    summarize(N = n(), Time.Sum = sum(Time.Sum), .groups = "drop_last") |>
    ungroup(Operation) |>
    arrange(Operation) |>
    mutate(TaskHeight = Time.Sum / (DurationTS * N)) |>
    mutate(TaskPosition = cumsum(TaskHeight) - TaskHeight) |>
    ggplot() +
    coord_cartesian(ylim = c(0, 1)) +
    geom_rect(alpha = 0.5, aes(
      fill = Operation,
      xmin = StartTS, xmax = EndTS,
      ymin = TaskPosition, ymax = TaskPosition + TaskHeight
    )) +
    ylab("Load [%]") +
    xlab("Time [seconds]") +
    load_theme +
    theme(legend.position = "none")
  p_rank / p_load + plot_layout(heights = c(1, 0.15))
}

generate_combined_load_chart <- function(real_df, sim_df) {
  load_theme <- list(
    theme_bw(base_size = 14),
    theme(
      panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_text(size = 10)
    )
  )
  facet <- facet_wrap(~src, labeller = as_labeller(c("real" = "Real", "sim" = "Simulation")))
  combined <- bind_rows(
    list(real = compute_load_data(real_df), sim = compute_load_data(sim_df)),
    .id = "src"
  ) |> mutate(src = as.factor(src))
  yconfm <- combined |>
    distinct(Rank) |>
    arrange(Rank) |>
    mutate(Position = 0:(n() - 1))
  p_rank <- combined |>
    filter(Time.Sum != 0) |>
    group_by(src, Rank, Slice, StartTS, EndTS, DurationTS, Operation) |>
    summarize(Time.Sum = sum(Time.Sum), .groups = "drop_last") |>
    ungroup(Operation) |>
    arrange(Operation) |>
    mutate(TaskHeight = Time.Sum / DurationTS) |>
    mutate(TaskPosition = cumsum(TaskHeight) - TaskHeight) |>
    left_join(yconfm, by = "Rank") |>
    ggplot() +
    scale_y_continuous(breaks = yconfm$Position, labels = yconfm$Rank) +
    geom_rect(alpha = 0.5, aes(
      fill = Operation,
      xmin = StartTS, xmax = EndTS,
      ymin = Position + TaskPosition,
      ymax = Position + TaskPosition + TaskHeight
    )) +
    facet +
    ylab("Rank [index]") +
    load_theme +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p_load <- combined |>
    filter(Time.Sum != 0) |>
    group_by(src, Slice, StartTS, EndTS, DurationTS, Operation) |>
    summarize(N = n(), Time.Sum = sum(Time.Sum), .groups = "drop_last") |>
    ungroup(Operation) |>
    arrange(Operation) |>
    mutate(TaskHeight = Time.Sum / (DurationTS * N)) |>
    mutate(TaskPosition = cumsum(TaskHeight) - TaskHeight) |>
    ggplot() +
    coord_cartesian(ylim = c(0, 1)) +
    geom_rect(alpha = 0.5, aes(
      fill = Operation,
      xmin = StartTS, xmax = EndTS,
      ymin = TaskPosition, ymax = TaskPosition + TaskHeight
    )) +
    facet +
    ylab("Load [%]") +
    xlab("Time [seconds]") +
    load_theme +
    theme(
      legend.position = "none", strip.text = element_blank(),
      strip.background = element_blank()
    )
  p_rank / p_load + plot_layout(heights = c(1, 0.15))
}

simplify_dataset <- function(df) {
  df <- df |>
    filter(problem_size %in% problem_sizes) |>
    arrange(Rank, problem_size, Start) |>
    group_by(Rank, problem_size) |>
    mutate(End = coalesce(min(End, lead(Start)), End)) |>
    mutate(Duration = End - Start) |>
    ungroup()
  core_start <- df |>
    filter(Operation == "MPI_Irecv") |>
    pull(Start) |>
    min()
  core_end <- df |>
    filter(Operation == "MPI_Waitall") |>
    pull(End) |>
    max()
  df <- df |>
    filter(Start < core_end) |>
    filter(End > core_start) |>
    mutate(Start = Start - core_start) |>
    mutate(End = End - core_start) |>
    filter(Operation %in% c("MPI_Irecv", "MPI_Isend", "MPI_Waitall")) |>
    mutate(Operation = as.factor(Operation))
  compute_rows <- df |>
    group_by(Rank, problem_size) |>
    mutate(next.start = lead(Start)) |>
    filter(!is.na(next.start)) |>
    transmute(Start = End, End = next.start, Duration = End - Start, Operation = "Compute") |>
    filter(Duration > 0) |>
    ungroup()
  bind_rows(df, compute_rows) |>
    arrange(Rank, problem_size, Start)
}

summarized_dataset <- function(df) {
  df |>
    filter(Operation != "Compute") |>
    group_by(Rank, problem_size, Operation) |>
    arrange(Start, .by_group = TRUE) |>
    mutate(event_id = row_number()) |>
    ungroup() |>
    group_by(problem_size, Operation, event_id) |>
    summarise(
      Duration = median(Duration),
      End = median(End),
      Start = median(Start),
      .groups = "drop"
    ) |>
    group_by(problem_size) |>
    summarize(mpi_time = sum(Duration), total_time = max(End) - min(Start)) |>
    mutate(compute_time = total_time - mpi_time) |>
    relocate(problem_size, compute_time, mpi_time, total_time)
}

filter_time_max <- function(df) {
  if (is.null(time_max)) {
    df
  } else {
    df |>
      filter(Start < time_max) |>
      mutate(End = pmin(time_max, End))
  }
}

save_load_charts <- function(real_dataset = NULL, simulation_dataset = NULL) {
  walk(problem_sizes, \(size) {
    base_dir <- "plots"
    dir.create(base_dir, showWarnings = FALSE)
    plot <- if (!is.null(real_dataset) && !is.null(simulation_dataset)) {
      generate_combined_load_chart(
        real_dataset |> filter(problem_size == size) |> filter_time_max(),
        simulation_dataset |> filter(problem_size == size) |> filter_time_max()
      )
    } else {
      df <- if (!is.null(real_dataset)) real_dataset else simulation_dataset
      generate_load_chart(compute_load_data(df |> filter(problem_size == size) |> filter_time_max()))
    }
    pdf_path <- paste0(base_dir, "/load_", size, ".pdf")
    ggsave(pdf_path, plot, width = 15, height = 8)
    print(paste0("Saving ", pdf_path, "..."))
  })
}

save_gantt_charts <- function(real_dataset = NULL, simulation_dataset = NULL) {
  walk(problem_sizes, \(size) {
    base_dir <- "plots"
    dir.create(base_dir, showWarnings = FALSE)
    if (!is.null(real_dataset) && !is.null(simulation_dataset)) {
      plot <- generate_combined_gantt_chart(
        real_dataset |> filter(problem_size == size) |> filter_time_max(),
        simulation_dataset |> filter(problem_size == size) |> filter_time_max()
      )
    } else {
      df <- if (!is.null(real_dataset)) real_dataset else simulation_dataset
      plot <- generate_gantt_chart(df |> filter(problem_size == size) |> filter_time_max())
    }
    pdf_path <- paste0(base_dir, "/gantt_", size, ".pdf")
    ggsave(pdf_path, plot, width = 15, height = 6)
    print(paste0("Saving ", pdf_path, "..."))
  })
}

real_dataset <- if (has_real) {
  load_real_traces(real_traces_path) |>
    filter(run_number == 1) |>
    select(-run_number) |>
    simplify_dataset()
} else {
  NULL
}

simulation_dataset <- if (has_sim) {
  load_simulation(sim_traces_path) |>
    simplify_dataset()
} else {
  NULL
}

save_gantt_charts(real_dataset, simulation_dataset)
save_load_charts(real_dataset, simulation_dataset)

if (has_real && has_sim) {
  real_summary <- real_dataset |> summarized_dataset()
  simulation_summary <- simulation_dataset |> summarized_dataset()
  comparison_summary <- inner_join(real_summary, simulation_summary, by = "problem_size") |>
    rename(
      real_compute = compute_time.x, real_mpi = mpi_time.x, real_total = total_time.x,
      sim_compute  = compute_time.y, sim_mpi  = mpi_time.y, sim_total  = total_time.y
    )
  print(comparison_summary)
}
