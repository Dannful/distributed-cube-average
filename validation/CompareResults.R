library(here)
library(digest)

ground_truth_path <- here::here("validation/ground_truth.dc")
predicted_path <- here::here("validation/predicted.dc")

read_floats <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  con <- file(file_path, "rb")
  on.exit(close(con))

  file_size <- file.info(file_path)$size
  num_floats <- file_size / 4

  floats <- readBin(con, what = "numeric", n = num_floats, size = 4)

  floats
}

# Read the float data from both files
ground_truth_floats <- read_floats(ground_truth_path)
predicted_floats <- read_floats(predicted_path)

# Compare the contents
if (length(ground_truth_floats) != length(predicted_floats)) {
  output <- "✖️ Validation error: the files have different sizes."
} else {
  # Compute differences
  differences <- ground_truth_floats - predicted_floats
  abs_differences <- abs(differences)

  if (all(abs_differences == 0)) {
    output <- "✅ Validation successful! The files are numerically identical."
  } else {
    max_abs_diff <- max(abs_differences)
    mean_abs_diff <- mean(abs_differences)
    num_diffs <- sum(abs_differences > 1e-9) # Count significant differences

    output <- paste0(
      "✖️ Validation error: the outputs differ numerically.\n",
      "Number of significant differences: ", num_diffs, "\n",
      "Maximum absolute difference: ", format(max_abs_diff, scientific = TRUE, digits = 4), "\n",
      "Mean absolute difference: ", format(mean_abs_diff, scientific = TRUE, digits = 4)
    )
  }
}

# Print the result
cat(output, "\n")
