library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Invalid usage. Use: <real traces dir> <simulation traces dir> <problem sizes>.")
}

# testing
args[1] <- "./analysis/2026-03-26/"
args[2] <- "./sim_traces/"
args[3] <- "128"

real_traces_path <- args[1]
sim_traces_path <- args[2]
problem_sizes <- as.numeric(str_split_1(args[3], ","))

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
    )
}

generate_combined_gantt_chart <- function(real_df, sim_df) {
  bind_rows(list(real = real_df, sim = sim_df), .id = "src") |>
    mutate(src = as.factor(src)) |>
    filter(Start < 3) |>
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
  df |>
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
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_text(size = 10)
    )
}

generate_combined_load_chart <- function(real_df, sim_df) {
  bind_rows(
    list(real = compute_load_data(real_df), sim = compute_load_data(sim_df)),
    .id = "src"
  ) |>
    mutate(src = as.factor(src)) |>
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
    facet_wrap(~src, labeller = as_labeller(c("real" = "Real", "sim" = "Simulation"))) +
    ylab("Load [%]") +
    xlab("Time [seconds]") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_text(size = 10)
    )
}

generate_rank_chart <- function(df) {
  yconfm <- df |>
    distinct(Rank) |>
    arrange(Rank) |>
    mutate(Position = 0:(n() - 1))
  df |>
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
    xlab("Time [seconds]") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_text(size = 10)
    )
}

generate_combined_rank_chart <- function(real_df, sim_df) {
  combined <- bind_rows(
    list(real = compute_load_data(real_df), sim = compute_load_data(sim_df)),
    .id = "src"
  ) |> mutate(src = as.factor(src))
  yconfm <- combined |>
    distinct(Rank) |>
    arrange(Rank) |>
    mutate(Position = 0:(n() - 1))
  combined |>
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
    facet_wrap(~src, labeller = as_labeller(c("real" = "Real", "sim" = "Simulation"))) +
    ylab("Rank [index]") +
    xlab("Time [seconds]") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_text(size = 10)
    )
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

save_load_charts <- function(real_dataset, simulation_dataset) {
  walk(problem_sizes, \(size) {
    base_dir <- "plots"
    dir.create(base_dir, showWarnings = FALSE)
    real_filtered <- real_dataset |> filter(problem_size == size)
    sim_filtered <- simulation_dataset |> filter(problem_size == size)
    p_rank <- generate_combined_rank_chart(real_filtered, sim_filtered) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    p_load <- generate_combined_load_chart(real_filtered, sim_filtered) +
      theme(legend.position = "none", strip.text = element_blank(), strip.background = element_blank())
    plot <- p_rank / p_load + plot_layout(heights = c(1, 0.15))
    pdf_path <- paste0(base_dir, "/load_", size, ".pdf")
    ggsave(pdf_path, plot, width = 15, height = 8)
    print(paste0("Saving ", pdf_path, "..."))
  })
}

save_gantt_charts <- function(real_dataset, simulation_dataset) {
  walk(problem_sizes, \(size) {
    base_dir <- "plots"
    dir.create(base_dir, showWarnings = FALSE)
    plot <- generate_combined_gantt_chart(real_dataset |>
      filter(problem_size == size), simulation_dataset |>
      filter(problem_size == size))
    pdf_path <- paste0(base_dir, "/gantt_", size, ".pdf")
    ggsave(pdf_path, plot, width = 15, height = 6)
    print(paste0("Saving ", pdf_path, "..."))
  })
}

real <- load_real_traces(real_traces_path) |>
  filter(run_number == 1) |>
  select(-run_number)

real_dataset <- load_real_traces(real_traces_path) |>
  filter(run_number == 1) |>
  select(-run_number) |>
  simplify_dataset()

simulation_dataset <- load_simulation(sim_traces_path) |>
  simplify_dataset()

save_gantt_charts(real_dataset, simulation_dataset)
save_load_charts(real_dataset, simulation_dataset)

real_summary <- real_dataset |>
  summarized_dataset()
simulation_summary <- simulation_dataset |>
  summarized_dataset()

comparison_summary <- inner_join(real_summary, simulation_summary, by = "problem_size") |>
  rename(
    real_compute = compute_time.x, real_mpi = mpi_time.x, real_total = total_time.x,
    sim_compute = compute_time.y, sim_mpi = mpi_time.y, sim_total = total_time.y
  )
print(comparison_summary)
