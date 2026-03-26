suppressMessages(library(tidyverse))
library(stringr)
library(grid)

# Setup directories
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("At least one argument must be supplied.", call. = FALSE)
}
base_dir <- args[1]
sim_traces_dir <- args[2]

# Helper to load and clean trace data into a dataframe
load_trace <- function(csv_path) {
    if (!file.exists(csv_path))
        return(NULL)
    df <- read_csv(csv_path, show_col_types = FALSE, col_names = c("Nature", "Container",
        "Type", "Start", "End", "Duration", "Imbrication", "Value")) |>
        mutate(Container = as.integer(gsub("rank-?", "", Container))) |>
        mutate(Value = str_trim(gsub("^\\s*PMPI_", "MPI_", str_trim(Value)))) |>
        mutate(Rank = Container) |>
        mutate(Operation = as.factor(Value)) |>
        filter(Operation %in% c("MPI_Irecv", "MPI_Isend", "MPI_Waitall"))
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
        mutate(Start = pmax(Start, core_start), End = pmin(End, core_end), Duration = End -
            Start)
    df |>
        mutate(Start = Start - min(df$Start)) |>
        mutate(End = End - min(df$Start))
}

# Function to create a ggplot Gantt chart
create_gantt_plot <- function(df_states, plot_title, x_lim = NULL) {
    if (is.null(df_states) || nrow(df_states) == 0)
        return(NULL)

    min_rank <- min(df_states$Rank, na.rm = TRUE)
    max_rank <- max(df_states$Rank, na.rm = TRUE)
    mid_breaks <- seq(min_rank + 0.5, max_rank + 0.5, by = 1)
    makespan <- max(df_states$End, na.rm = TRUE)
    epsilon <- makespan/200
    p <- df_states |>
        filter(!str_detect(Operation, "MPI_Init|MPI_Finalize|MPI_Allgather|MPI_Comm|MPI_Cart|MPI_Dims")) |>
        ggplot() + geom_hline(yintercept = seq(min_rank + 0.5, max_rank + 1 - 0.5,
        by = 1), color = "grey85", linewidth = 0.1) + geom_rect(aes(xmin = Start,
        xmax = pmax(End, Start + epsilon), ymin = Rank, ymax = Rank + 0.9, fill = Operation)) +
        scale_y_continuous(breaks = mid_breaks, labels = seq(min_rank, max_rank,
            by = 1)) + xlab("Time [seconds]") + ylab("Rank [count]") + scale_fill_manual(values = c(MPI_Irecv = "#FCA5A7",
        MPI_Isend = "#A6C8E5", MPI_Recv = "#E41A1C", MPI_Send = "#377EB8", MPI_Waitall = "#FF7F00",
        MPI_Barrier = "#4DAF4A"), breaks = c("MPI_Recv", "MPI_Irecv", "MPI_Send",
        "MPI_Isend", "MPI_Waitall", "MPI_Barrier")) + geom_text(data = data.frame(x = makespan,
        y = (max_rank - min_rank)/2), aes(x = x, y = y, label = paste0(round(makespan,
        3))), inherit.aes = FALSE, angle = 90, size = 2.5) + theme_bw(base_size = 10) +
        guides(fill = guide_legend(nrow = 1)) + theme(plot.margin = unit(c(0.1, 0.1,
        0.1, 0.1), "cm"), legend.margin = margin(t = 0, unit = "cm"), panel.grid = element_blank(),
        legend.position = "top", legend.justification = "left", legend.box.spacing = unit(1,
            "pt"), legend.box.margin = margin(0, 0, 0, 0), legend.title = element_blank())
    if (!is.null(x_lim)) {
        p <- p + coord_cartesian(xlim = x_lim)
    }
    p
}

# 1. Process Real Data
print("Loading real traces...")
rst_files <- list.files(path = base_dir, pattern = "\\.rst$", recursive = TRUE, full.names = TRUE)
run_dirs <- unique(dirname(rst_files))

df_real_all <- map_dfr(run_dirs, function(dir_path) {
    matches <- str_match(dir_path, "analysis/[^/]+/(?<run>[0-9]+)/(?<size>[0-9]+)")
    run_num <- as.integer(matches[1, "run"])
    size_num <- as.integer(matches[1, "size"])

    trace_path <- tempfile(fileext = ".trace")
    csv_path <- tempfile(fileext = ".csv")
    system(paste0("aky_converter -l ", dir_path, "/*.rst > ", trace_path, " 2>/dev/null"))
    system(paste0("pj_dump -z -l 9 ", trace_path, " | grep ^State > ", csv_path,
        " 2>/dev/null"))

    df <- load_trace(csv_path) |>
        mutate(run_number = run_num, problem_size = size_num)
    unlink(c(trace_path, csv_path))
    df
})

df_real_summary <- df_real_all |>
    group_by(problem_size, run_number, Rank) |>
    summarise(total_time = max(End, na.rm = TRUE) - min(Start, na.rm = TRUE), mpi_time = sum(Duration,
        na.rm = TRUE), .groups = "drop") |>
    mutate(compute_time = total_time - mpi_time) |>
    group_by(problem_size) |>
    summarise(total_time = mean(total_time), mpi_time = mean(mpi_time), compute_time = mean(compute_time))

# 2. Process Simulated Data Load all simulation trace data to calculate
# MPI/compute time breakdown
print("Loading simulation traces...")
size_dirs <- list.dirs(sim_traces_dir, full.names = FALSE, recursive = FALSE)
df_sim_all_traces <- map_dfr(size_dirs, function(sz) {
    sim_csv <- file.path(sim_traces_dir, sz, "dc.csv")
    if (file.exists(sim_csv)) {
        load_trace(sim_csv) |>
            mutate(problem_size = as.integer(sz))
    } else {
        NULL
    }
})

df_sim_summary <- df_sim_all_traces |>
    group_by(problem_size, Rank) |>
    summarise(total_time = max(End, na.rm = TRUE) - min(Start, na.rm = TRUE), mpi_time = sum(Duration,
        na.rm = TRUE), .groups = "drop") |>
    mutate(compute_time = total_time - mpi_time) |>
    group_by(problem_size) |>
    summarise(total_time = mean(total_time), mpi_time = mean(mpi_time), compute_time = mean(compute_time))

if (nrow(df_sim_summary) > 0 && nrow(df_real_summary) > 0) {
    # 3. Compare and Plot
    merged_summary <- inner_join(df_real_summary, df_sim_summary, by = "problem_size",
        suffix = c("_real", "_sim"))
    if (nrow(merged_summary) > 0) {
        dir.create("plots", showWarnings = FALSE)

        # Print comparison table to stdout
        cat("\nTime Breakdown Comparison (seconds):\n")
        comparison_df <- merged_summary |>
            select(problem_size, real_compute = compute_time_real, real_mpi = mpi_time_real,
                sim_compute = compute_time_sim, sim_mpi = mpi_time_sim)
        print(as.data.frame(comparison_df))
        cat("\n")

        # Prepare data for stacked bar chart
        plot_data <- merged_summary |>
            select(problem_size, compute_time_real, mpi_time_real, compute_time_sim,
                mpi_time_sim) |>
            pivot_longer(cols = -problem_size, names_to = c("time_type", "source"),
                names_pattern = "(.+)_(real|sim)", values_to = "time") |>
            mutate(source = factor(ifelse(source == "real", "Real", "Simulated"),
                levels = c("Real", "Simulated")), time_type = factor(ifelse(time_type ==
                "compute_time", "Computation", "MPI"), levels = c("MPI", "Computation")))

        # Stacked bar comparison plot
        p_comp <- ggplot(plot_data, aes(x = factor(problem_size), y = time, fill = time_type)) +
            geom_bar(stat = "identity", position = "stack") + facet_wrap(~source) +
            scale_fill_manual(values = c(Computation = "#4DAF4A", MPI = "#E41A1C")) +
            labs(title = "Execution Time Breakdown: Computation vs MPI", x = "Problem Size",
                y = "Time (s)", fill = "Time Type") + theme_minimal() + theme(legend.position = "top")
        ggsave("plots/comparison_plot.pdf", p_comp, width = 10, height = 6)

        # Gantt Chart Comparisons
        for (sz in merged_summary$problem_size) {
            # Compute mean real Gantt data
            df_r_mean <- df_real_all |>
                filter(problem_size == sz) |>
                group_by(run_number, Rank) |>
                arrange(Start) |>
                mutate(event_id = row_number()) |>
                group_by(Rank, event_id, Value) |>
                summarise(Start = mean(Start), End = mean(End), Duration = mean(Duration),
                  .groups = "drop") |>
                mutate(Operation = fct_drop(as.factor(Value)))

            # Filter sim traces for this size (already loaded above)
            df_sim_size <- df_sim_all_traces |>
                filter(problem_size == sz)

            if (!is.null(df_sim_size) && nrow(df_sim_size) > 0) {
                # Compute mean sim Gantt data
                df_s_mean <- df_sim_size |>
                  group_by(Rank) |>
                  arrange(Start) |>
                  mutate(event_id = row_number()) |>
                  group_by(Rank, event_id, Value) |>
                  summarise(Start = mean(Start), End = mean(End), Duration = mean(Duration),
                    .groups = "drop") |>
                  mutate(Operation = fct_drop(as.factor(Value)))

                max_time <- max(c(df_r_mean$End, df_s_mean$End), na.rm = TRUE)
                lims <- c(0, max_time)

                p_r <- create_gantt_plot(df_r_mean, paste0("Real Data (Mean) - Size ",
                  sz), x_lim = lims)
                p_s <- create_gantt_plot(df_s_mean, paste0("Simulated Data (Mean) - Size ",
                  sz), x_lim = lims)

                if (!is.null(p_r) && !is.null(p_s)) {
                  pdf(file.path("plots", paste0("combined_gantt_size_", sz, ".pdf")),
                    width = 10, height = 5)
                  grid.newpage()
                  pushViewport(viewport(layout = grid.layout(2, 1)))
                  print(p_s, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
                  print(p_r, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
                  dev.off()
                }
            }
        }
        print("Plots generated in plots/ directory.")
    }
}

