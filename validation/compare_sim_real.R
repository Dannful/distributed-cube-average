suppressMessages(library(tidyverse))
library(stringr)

base_dir <- "analysis/2026-01-20"
if (!dir.exists(base_dir)) {
    warning(paste("Real data directory not found:", base_dir, "- Skipping real data loading."))
    df_real <- data.frame()
} else {
    files <- list.files(path = base_dir, pattern = "dc\.mpiP$", recursive = TRUE, full.names = TRUE)
    
    parse_mpip <- function(file_path) {
        path_parts <- str_match(file_path, "/([0-9]+)/([0-9]+)/dc\.mpiP$")
        if (any(is.na(path_parts))) return(NULL)
        
        run_number <- as.integer(path_parts[2])
        problem_size <- as.integer(path_parts[3])
        lines <- readLines(file_path, warn = FALSE)
        time_idx <- grep("@.*Time.*seconds", lines)
        
        if (length(time_idx) > 0) {
            subset_lines <- lines[time_idx:min(time_idx + 50, length(lines))]
            agg_line_idx <- grep("^\s*\*\s+[0-9.]+\s+[0-9.]+\s+([0-9.]+)$", subset_lines)
            parsed_data <- NULL
            
            if (length(agg_line_idx) > 0) {
                agg_line <- subset_lines[agg_line_idx[1]]
                parsed_data <- str_match(agg_line, "^\s*\*\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)$")
            } else {
                task0_line_idx <- grep("^\s*0\s+[0-9.]+\s+[0-9.]+\s+([0-9.]+)$", subset_lines)
                if (length(task0_line_idx) > 0) {
                    agg_line <- subset_lines[task0_line_idx[1]]
                    parsed_data <- str_match(agg_line, "^\s*0\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)$")
                }
            }
            
            if (!is.null(parsed_data) && !any(is.na(parsed_data))) {
                total_time <- as.numeric(parsed_data[2])
                mpi_time <- as.numeric(parsed_data[3])
                computation_time <- total_time - mpi_time
                return(data.frame(run = run_number, problem_size = problem_size, 
                                  mpi_time = mpi_time, computation_time = computation_time, 
                                  total_time = total_time, type = "Real"))
            }
        }
        NULL
    }
    
    df_real <- map_dfr(files, parse_mpip)
}

args <- commandArgs(trailingOnly = TRUE)
sim_file <- if (length(args) > 0) args[1] else "simulation_results.csv"

if (!file.exists(sim_file)) {
    stop(paste("Simulation results file not found:", sim_file))
}

df_sim <- read_csv(sim_file, show_col_types = FALSE) |>
    mutate(type = "Simulation")

cols <- c("run", "problem_size", "mpi_time", "computation_time", "total_time", "type")
df_combined <- bind_rows(
    df_real |> select(all_of(cols)),
    df_sim |> select(all_of(cols))
)

df_summary <- df_combined |>
    group_by(type, problem_size) |>
    summarise(
        mpi_time = mean(mpi_time), 
        computation_time = mean(computation_time),
        total_time = mean(total_time),
        .groups = "drop"
    ) |>
    mutate(mpi_percentage = mpi_time / total_time)

p <- ggplot(data = df_summary, mapping = aes(x = problem_size, y = mpi_percentage, color = type)) + 
    geom_line() + 
    geom_point() + 
    scale_x_continuous(trans='log2') +
    labs(title = "MPI Percentage: Real vs Simulation", 
         x = "Problem Size (N)", y = "MPI Time / Total Time") +
    theme_minimal()

ggsave("comparison_chart.pdf", p, width = 8, height = 6)
print("Comparison chart saved to comparison_chart.pdf")

df_wide <- df_summary |>
    select(problem_size, type, mpi_percentage) |>
    pivot_wider(names_from = type, values_from = mpi_percentage) |>
    na.omit()

if (nrow(df_wide) > 2) {
    model <- lm(Real ~ Simulation, data = df_wide)
    r_squared <- summary(model)$r.squared
    
    print(paste("R-squared (MPI Percentage):", round(r_squared, 4)))
    
    df_wide_time <- df_summary |>
        select(problem_size, type, total_time) |>
        pivot_wider(names_from = type, values_from = total_time) |>
        na.omit()
    
    model_time <- lm(Real ~ Simulation, data = df_wide_time)
    r_squared_time <- summary(model_time)$r.squared
    print(paste("R-squared (Total Time):", round(r_squared_time, 4)))
    
} else {
    print("Not enough data points to calculate R-squared.")
}
