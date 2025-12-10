# Read the data
suppressMessages(library(tidyverse))
read_csv("dc.csv",
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
  # Create the nice MPI rank and operations identifiers
  mutate(Container = as.integer(gsub("rank-", "", Container))) |>
  mutate(Value = gsub("^PMPI_", "MPI_", Value)) |>
  # Rename some columns so it can better fit MPI terminology
  mutate(Rank = Container) |>
  mutate(Operation = as.factor(Value)) -> df.states

# Let's disregard this for now

# df.states |>
#    filter(Operation == "MPI_Waitall") |>
#    filter(End == max(End)) |>
#    pull(End) -> df.end

# df.states |>
#    filter(Operation == "MPI_Send") |>
#    filter(End < df.end) |>
#    filter(End == max(End)) |>
#    pull(End) -> df.start

# df.states <- df.states |>
#    filter(Start > df.start) |>
#    filter(End < df.end)

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

df.overall <- df.states |>
  filter(Rank == 0, Duration > 0)
true_start <- max(df.overall |>
  filter(Operation == "MPI_Send") |>
  select(Start))
df.overall <- df.overall |>
  filter(Start >= true_start) |>
  select(Start, End, Duration)
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

# Draw the Gantt Chart
df.states |>
  filter(Operation != "MPI_Allgather") |> # does not appear in the plot
  filter(Operation != "MPI_Finalize") |> # does not appear in the plot
  filter(Operation != "MPI_Init") |> # does not appear in the plot
  ggplot() + # Each MPI operation is becoming a rectangle ggplot() +
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
    ) # ,
    # color = "grey90",     # thin subtle border
    # linewidth = 0.1        # very thin line
  ) +
  scale_y_continuous(
    breaks = mid_breaks,
    labels = seq(min_rank, max_rank, by = 1)
  ) +
  xlab("Time [seconds]") +
  ylab("Rank [count]") +
  # scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(
    values = c(
      "MPI_Irecv" = "#FCA5A7", # light red
      "MPI_Isend" = "#A6C8E5", # light blue
      "MPI_Recv" = "#E41A1C", # red
      "MPI_Send" = "#377EB8", # blue
      "MPI_Waitall" = "#FF7F00" # orange
    ),
    breaks = c(
      "MPI_Recv",
      "MPI_Irecv",
      "MPI_Send",
      "MPI_Isend",
      "MPI_Waitall"
    )
  ) +
  # vertical makespan label (drawn once)
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

# Save the plot in a PDF file (dimensions in inches)
ggsave("smpi.pdf", plot, width = 10, height = 2)

# Save the plot in a PNG file
ggsave("smpi.png", plot, width = 10, height = 2, dpi = 300)

library(plotly)
interactive_plot <- ggplotly(plot)

library(htmlwidgets)

saveWidget(
  widget = interactive_plot,
  file = "smpi.html",
  selfcontained = TRUE
)
