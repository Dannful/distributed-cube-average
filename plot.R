# Read the data
suppressMessages(library(tidyverse))
read_csv("dc.csv",
         show_col_types = FALSE,
         col_names = c("Nature",
                       "Container",
                       "Type",
                       "Start",
                       "End",
                       "Duration",
                       "Imbrication",
                       "Value")) |>
    select(-Type) |>
    select(-Imbrication) |>
    # Create the nice MPI rank and operations identifiers
    mutate(Container = as.integer(gsub("rank-", "", Container))) |>
    mutate(Value = gsub("^PMPI_", "MPI_", Value)) |>
    # Rename some columns so it can better fit MPI terminology
    mutate(Rank = Container) |>
    mutate(Operation = as.factor(Value)) -> df.states

# Let's disregard this for now

#df.states |>
#    filter(Operation == "MPI_Waitall") |>
#    filter(End == max(End)) |>
#    pull(End) -> df.end

#df.states |>
#    filter(Operation == "MPI_Send") |>
#    filter(End < df.end) |>
#    filter(End == max(End)) |>
#    pull(End) -> df.start

# df.states <- df.states |>
#    filter(Start > df.start) |>
#    filter(End < df.end)

df.states |>
    select(-Nature, -Container, -Duration) |>
    group_by(Rank) |>
    arrange(Start) |>
    mutate(C.Start = End) |>
    mutate(C.End = lead(Start)) |>
    mutate(C.Duration = C.End - C.Start) |>
    summarize(Compute.Cost = sum(C.Duration, na.rm=TRUE))

# Draw the Gantt Chart
df.states |>
    ggplot() + # Each MPI operation is becoming a rectangle ggplot() +
    geom_rect(aes(xmin = Start, xmax = End, ymin = Rank, ymax = Rank + 0.9, fill = Operation)) +
    # Cosmetics
    xlab("Time [seconds]") +
    ylab("Rank [count]") +
    theme_bw(base_size = 14) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.margin = margin(t = 0, unit = "cm"),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.justification = "left",
          legend.box.spacing = unit(0, "pt"),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 10)) -> plot

# Save the plot in a PDF file (dimensions in inches)
ggsave("smpi.pdf", plot, width = 10, height = 3)
