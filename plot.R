# Read the data
suppressMessages(library(tidyverse))
library(plotly)
library(htmlwidgets)

# Argument handling
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "dc.csv"

if (!file.exists(input_file)) {
  stop(paste("File not found:", input_file))
}

# Determine file type
# Read file to detect format. 
# We read the whole file (or enough of it) because mpiP output might be at the end of a log.
lines <- readLines(input_file)
is_mpip <- any(grepl("^@ MPIP Summary", lines)) || any(grepl("^@ Command :", lines))

if (is_mpip) {
  # --- mpiP Parsing ---
  print(paste("Detected mpiP output in:", input_file))
  
  # Find Time section
  # @ Time (seconds)
  # @   App    MPI    %
  time_idx <- grep("^@ Time \\(seconds\\)", lines)
  
  if (length(time_idx) > 0) {
    # Extract lines that look like data: "@   1.23   0.45  36.59"
    # We search in the lines following the header
    subset_lines <- lines[time_idx:length(lines)]
    data_lines <- subset_lines[grep("^@ \\s*[0-9.]+\\s+[0-9.]+\\s+[0-9.]+", subset_lines)]
    
    if (length(data_lines) > 0) {
      # Parse the first data line found (assuming aggregate or representative)
      parsed_data <- str_match(data_lines[1], "^@ \\s*([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)")
      
      if (!any(is.na(parsed_data))) {
        total_time <- as.numeric(parsed_data[2])
        mpi_time <- as.numeric(parsed_data[3])
        
        computation_time <- total_time - mpi_time
        computation_percentage <- computation_time / total_time
        
        print(paste0("Total time: ", total_time))
        print(paste0("MPI time: ", mpi_time))
        print(paste0("Computation time: ", computation_time))
        print(paste0("Total coverage: ", computation_percentage))
      } else {
        print("Error parsing data line in Time section.")
      }
    } else {
      print("No data lines found in Time section.")
    }
  } else {
    print("No Time section found in mpiP output.")
  }

} else {
  # --- Paje/CSV Parsing ---
  print(paste("Detected CSV/Paje output in:", input_file))
  
  read_csv(input_file,
    show_col_types = FALSE,
    col_names = c(
      "Nature",
      "Container",
      "Type",
      "Start",
      "End",
      "Duration",
      "Imbrication",
      "Value"
    )
  ) |>
    select(-Type) |>
    select(-Imbrication) |>
    mutate(Container = as.integer(gsub("rank-", "", Container))) |>
    mutate(Value = gsub("^PMPI_", "MPI_", Value)) |>
    mutate(Rank = Container) |>
    mutate(Operation = as.factor(Value)) -> df.states

  # The "Compute Cost" per rank
  df.states |>
    select(-Nature, -Container, -Duration) |>
    group_by(Rank) |>
    arrange(Start) |>
    mutate(C.Start = End) |>
    mutate(C.End = lead(Start)) |>
    mutate(C.Duration = C.End - C.Start) |>
    summarize(Compute.Cost = sum(C.Duration, na.rm = TRUE)) |>
    pull(Compute.Cost) |>
    summary()

  min_rank <- min(df.states$Rank)
  max_rank <- max(df.states$Rank)
  mid_breaks <- seq(min_rank + 0.5, max_rank + 0.5, by = 1)
  makespan <- max(df.states$End)

  df.overall <- df.states |
    filter(Rank == 0, Duration > 0)
  
  if (nrow(df.overall) > 0) {
    true_start <- max(df.overall |
      filter(Operation == "MPI_Send") |
      select(Start))
      
    if (is.infinite(true_start) || is.na(true_start) || nrow(df.overall |> filter(Operation == "MPI_Send")) == 0) {
        true_start <- min(df.overall$Start)
    }

    df.overall <- df.overall |
      filter(Start >= true_start) |>
      select(Start, End, Duration)
    
    if (nrow(df.overall) > 0) {
        min_time <- min(df.overall$Start)
        max_time <- max(df.overall$End)
        total_time <- max_time - min_time
        mpi_time <- sum(df.overall$Duration)
        computation_time <- total_time - mpi_time
        computation_percentage <- computation_time / total_time
        print(paste0("Total time: ", total_time))
        print(paste0("MPI time: ", mpi_time))
        print(paste0("Computation time: ", computation_time))
        print(paste0("Total coverage: ", computation_percentage))
    }
  }

  # Draw the Gantt Chart
  df.states |
    filter(Operation != "MPI_Allgather") |> 
    filter(Operation != "MPI_Finalize") |> 
    filter(Operation != "MPI_Init") |> 
    ggplot() + 
    geom_hline(
      yintercept = seq(min(df.states$Rank) + 0.5,
        max(df.states$Rank) + 1 - 0.5,
        by = 1
      ),
      color = "grey85",
      linewidth = 0.1
    ) + 
    geom_rect(
      aes(
        xmin = Start, xmax = End,
        ymin = Rank, ymax = Rank + 0.9,
        fill = Operation
      ) 
    ) + 
    scale_y_continuous(
      breaks = mid_breaks,
      labels = seq(min_rank, max_rank, by = 1)
    ) + 
    xlab("Time [seconds]") +
    ylab("Rank [count]") +
    scale_fill_manual(
      values = c(
        "MPI_Irecv" = "#FCA5A7", 
        "MPI_Isend" = "#A6C8E5", 
        "MPI_Recv" = "#E41A1C", 
        "MPI_Send" = "#377EB8", 
        "MPI_Waitall" = "#FF7F00" 
      ),
      breaks = c(
        "MPI_Recv",
        "MPI_Irecv",
        "MPI_Send",
        "MPI_Isend",
        "MPI_Waitall"
      )
    ) + 
    geom_text(
      data = data.frame(x = makespan, y = (max_rank - min_rank) / 2),
      aes(x = x, y = y, label = paste0(round(makespan, 0))),
      inherit.aes = FALSE,
      angle = 90,
      size = 2.5
    ) + 
    theme_bw(base_size = 10) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.margin = margin(t = 0, unit = "cm"),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.box.spacing = unit(1, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title = element_blank()
    ) -> plot

  ggsave("smpi.pdf", plot, width = 10, height = 2)
  ggsave("smpi.png", plot, width = 10, height = 2, dpi = 300)

  interactive_plot <- ggplotly(plot)

  saveWidget(
    widget = interactive_plot,
    file = "smpi.html",
    selfcontained = TRUE
  )
}
