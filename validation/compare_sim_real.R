suppressMessages(library(tidyverse))
library(stringr)

base_dir <- "analysis/2026-01-20"

if (!dir.exists(base_dir)) {
    stop(paste("Directory not found:", base_dir))
}

files <- list.files(path = base_dir, pattern = "dc\\.output$", recursive = TRUE, full.names = TRUE)

if (length(files) == 0) {
    warning("No dc.output files found.")
}

parse_output <- function(file_path) {
    path_parts <- str_match(file_path, "/([0-9]+)/([0-9]+)/dc\\.output$")

    if (any(is.na(path_parts))) {
        warning(paste("Could not parse path structure for:", file_path))
        return(NULL)
    }

    run_number <- as.integer(path_parts[2])
    problem_size <- as.integer(path_parts[3])

    lines <- readLines(file_path, warn = FALSE)
    header_idx <- grep("^rank,total_time,msamples_per_s", lines)

    if (length(header_idx) > 0) {
        csv_lines <- lines[header_idx[1]:length(lines)]
        df <- read.csv(text = paste(csv_lines, collapse = "\n"), colClasses=c("rank"="character"))
        df_rank0 <- df[df$rank == "*", ]
        if (nrow(df_rank0) > 0) {
            return(data.frame(run = run_number, problem_size = problem_size,
                total_time = df_rank0$total_time[1], msamples_per_s = df_rank0$msamples_per_s[1]))
        }
    }

    NULL
}

df_real <- map_dfr(files, parse_output)
if (nrow(df_real) > 0) {
    df_real <- df_real |>
        group_by(problem_size) |>
        summarise(total_time = mean(total_time))

    print(df_real)
} else {
    print("No real data parsed.")
}

sim_file <- here::here("simulation_results.csv")
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
