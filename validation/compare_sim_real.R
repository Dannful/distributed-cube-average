suppressMessages(library(tidyverse))
library(stringr)

base_dir <- "analysis/2026-01-20"

if (!dir.exists(base_dir)) {
    stop(paste("Directory not found:", base_dir))
}

files <- list.files(path = base_dir, pattern = "dc\\.mpiP$", recursive = TRUE, full.names = TRUE)

if (length(files) == 0) {
    stop("No dc.mpiP files found.")
}

parse_mpip <- function(file_path) {
    path_parts <- str_match(file_path, "/([0-9]+)/([0-9]+)/dc\\.mpiP$")

    if (any(is.na(path_parts))) {
        warning(paste("Could not parse path structure for:", file_path))
        return(NULL)
    }

    run_number <- as.integer(path_parts[2])
    problem_size <- as.integer(path_parts[3])

    lines <- readLines(file_path, warn = FALSE)

    time_idx <- grep("@.*Time.*seconds", lines)

    if (length(time_idx) > 0) {
        subset_lines <- lines[time_idx:min(time_idx + 50, length(lines))]

        agg_line_idx <- grep("^\\s*\\*\\s+[0-9.]+\\s+[0-9.]+\\s+([0-9.]+)$", subset_lines)

        parsed_data <- NULL

        if (length(agg_line_idx) > 0) {
            agg_line <- subset_lines[agg_line_idx[1]]
            parsed_data <- str_match(agg_line, "^\\s*\\*\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)$")
        } else {
            task0_line_idx <- grep("^\\s*0\\s+[0-9.]+\\s+[0-9.]+\\s+([0-9.]+)$",
                subset_lines)
            if (length(task0_line_idx) > 0) {
                agg_line <- subset_lines[task0_line_idx[1]]
                parsed_data <- str_match(agg_line, "^\\s*0\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)$")
            }
        }

        if (!is.null(parsed_data) && !any(is.na(parsed_data))) {
            total_time <- as.numeric(parsed_data[2])
            mpi_time <- as.numeric(parsed_data[3])
            computation_time <- total_time - mpi_time

            return(data.frame(run = run_number, problem_size = problem_size, mpi_time = mpi_time,
                computation_time = computation_time, total_time = total_time))
        }
    }

    NULL
}

df_real <- map_dfr(files, parse_mpip)
df_real <- df_real |>
    group_by(problem_size) |>
    summarise(mpi_time = mean(mpi_time), computation_time = mean(computation_time),
        total_time = mean(total_time)) |>
    mutate(mpi_percentage = mpi_time/total_time)

df_real

df_sim <- read.csv(here::here("simulation_results.csv"))
df_sim <- df_sim |>
    group_by(problem_size) |>
    summarise(mpi_time = mean(mpi_time), computation_time = mean(computation_time),
        total_time = mean(total_time)) |>
    mutate(mpi_percentage = mpi_time/total_time)
df_sim

model <- lm(df_real$mpi_percentage ~ df_sim$mpi_percentage)
summary(model)
