suppressMessages(library(tidyverse))
df <- read_csv("dc.state.csv", progress=FALSE, show_col_types=FALSE)
if (nrow(df) > 0) {
    total_time <- max(df$End) - min(df$Start)
    rank0 <- df[grep("rank-0", df$Container), ]
    print(paste0("Total time: ", total_time))
    print(paste0("MPI time: ", sum(rank0$Duration)))
    print(paste0("Computation time: ", total_time - sum(rank0$Duration)))
}
