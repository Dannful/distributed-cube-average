library(dplyr)
library(tibble)
library(optparse)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)

option_list <- list(
  make_option("--real-data-dir", type = "character", default = NULL,
    help = "Comma-separated list of real trace directories"),
  make_option("--simulation-data-dir", type = "character", default = NULL,
    help = "Comma-separated list of simulation trace directories")
)
opts <- parse_args(OptionParser(option_list = option_list))

parse_dir_args <- function(arg) {
  entries <- str_split_1(arg, ",")
  parts <- str_split_fixed(entries, ":", 2)
  set_names(parts[, 1], parts[, 2])
}
real_traces_paths <- if (!is.null(opts$`real-data-dir`)) parse_dir_args(opts$`real-data-dir`) else NULL
sim_traces_paths <- if (!is.null(opts$`simulation-data-dir`)) parse_dir_args(opts$`simulation-data-dir`) else NULL
has_real <- !is.null(real_traces_paths)
has_sim <- !is.null(sim_traces_paths)

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
    mutate(problem_size = as.integer(str_extract(source_file, "(?<=/)\\d+(?=/[^/]*$)"))) |>
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

simplify_dataset <- function(df) {
  df <- df |>
    arrange(Rank, problem_size, Start)
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
    filter(Operation %in% c("MPI_Irecv", "MPI_Isend", "MPI_Waitall"))
  compute_rows <- df |>
    group_by(Rank, problem_size) |>
    mutate(next.start = lead(Start)) |>
    filter(!is.na(next.start)) |>
    transmute(Start = End, End = next.start, Duration = End - Start, Operation = "Compute") |>
    filter(Duration > 0) |>
    ungroup()
  df <- bind_rows(df, compute_rows) |>
    arrange(Rank, problem_size, Start) |>
    mutate(Operation = as.factor(Operation))
  df |>
    group_by(problem_size) |>
    summarise(msamples = (first(problem_size)^3 * 1001) / ((max(End) - min(Start)) * 1e6))
}

base_fletcher <- tibble(
  problem_size = c(128, 256, 512),
  msamples = c(2635, 3773, 3870)
)

real_datasets <- if (has_real) {
  map_dfr(real_traces_paths, \(path) {
    load_real_traces(path) |>
      filter(run_number == 1) |>
      select(-run_number) |>
      simplify_dataset()
  }, .id = "source")
} else {
  NULL
}

simulation_datasets <- if (has_sim) {
  map_dfr(sim_traces_paths, \(path) {
    load_simulation(path) |>
      simplify_dataset()
  }, .id = "source")
} else {
  NULL
}

combined <- bind_rows(real_datasets, simulation_datasets, base_fletcher |> mutate(source = "Single-node"))

ggplot(combined, aes(x = problem_size, y = msamples, color = source)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = round(msamples)), vjust = -1, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Problem Size", y = "MSamples/s", color = "Source") +
  theme_minimal()

ggsave("msamples.pdf", width = 8, height = 5)

