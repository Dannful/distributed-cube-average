suppressMessages(library(tidyverse))
library(stringr)

base_dir <- "analysis/2026-03-06"

if (!dir.exists(base_dir)) {
    stop(paste("Directory not found:", base_dir))
}

files <- list.files(path = base_dir, pattern = "\\.rst$", recursive = TRUE, full.names = TRUE)

if (length(files) == 0) {
    warning("No dc.output files found.")
}

parse_output <- function(file_path) {
    matches <- str_match_all(file_path, "analysis/[^/]+/(?<run>[0-9]+)/(?<size>[0-9]+)")[[1]]
    run_number <- as.integer(matches[2])
    problem_size <- as.integer(matches[3])
    trace_path <- paste0(file_path, ".trace")
    csv_path <- paste0(file_path, ".csv")
    system(paste0("aky_converter -l ", file_path, " > ", trace_path))
    system(paste0("pj_dump -z -l 9 ", trace_path, " | grep ^State > ", csv_path))
    dataset <- read_csv(csv_path, show_col_types = FALSE, col_names = c("Nature",
        "Container", "Type", "Start", "End", "Duration", "Imbrication", "Value")) |>
        select(-Type) |>
        select(-Imbrication) |>
        mutate(Container = as.integer(gsub("rank-?", "", Container))) |>
        mutate(Value = gsub("^PMPI_", "MPI_", Value)) |>
        mutate(Rank = Container) |>
        mutate(Operation = as.factor(Value)) |>
        mutate(problem_size = problem_size) |>
        mutate(run_number = run_number) |>
        mutate(total_time = End - Start)
    system(paste0("rm ", trace_path))
    system(paste0("rm ", csv_path))
    dataset
}

print("Loading dataset...")
df_real <- map_dfr(files, parse_output)
print("Dataset loaded.")
if (nrow(df_real) > 0) {
    df_real <- df_real |>
        group_by(problem_size, run_number) |>
        summarise(total_time = max(End, na.rm = TRUE) - min(Start, na.rm = TRUE), .groups = "drop") |>
        group_by(problem_size) |>
        summarise(total_time = mean(total_time))

    print(df_real)
} else {
    print("No real data parsed.")
}

args <- commandArgs(trailingOnly = TRUE)
sim_file <- args[1]
if (file.exists(sim_file)) {
    df_sim <- read.csv(sim_file)
    df_sim <- df_sim |>
        group_by(problem_size) |>
        summarise(total_time = mean(total_time))
    print(df_sim)

    if (nrow(df_real) > 0 && nrow(df_sim) > 0) {
        model <- lm(df_real$total_time ~ df_sim$total_time)
        print(summary(model))
    }
} else {
    print("Simulation results file not found.")
}
