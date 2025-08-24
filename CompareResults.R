library(here)
library(digest)

ground_truth_path <- here::here("ground_truth.dc")
predicted_path <- here::here("predicted.dc")

ground_truth_hash <- digest::digest(file = ground_truth_path, algo = "sha256")
predicted_hash <- digest::digest(file = predicted_path, algo = "sha256")

output <- if (predicted_hash == ground_truth_hash) "✅ Validation successful!" else "✖️ Validation error: the outputs differ."
output
