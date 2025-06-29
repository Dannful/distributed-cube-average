#include "../include/log.h"
#include "unistd.h"
#include <stdarg.h>
#include <stdio.h>

void log_info(int rank, char *message,...) {
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  char worker_type[12] = {0};
  if(rank == 0) {
    snprintf(worker_type, sizeof(worker_type), "Coordinator");
  } else {
    snprintf(worker_type, sizeof(worker_type), "Worker %d", rank);
  }
  va_list args;
  va_start(args, message);
  printf("[INFO] %s - %s: ", hostname, worker_type);
  printf(message, args);
  printf("\n");
}

void log_error(int rank, char *message, ...) {
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  char worker_type[12] = {0};
  if(rank == 0) {
    snprintf(worker_type, sizeof(worker_type), "Coordinator");
  } else {
    snprintf(worker_type, sizeof(worker_type), "Worker %d", rank);
  }
  va_list args;
  va_start(args, message);
  fprintf(stderr, "[ERROR] %s - %s: ", hostname, worker_type);
  fprintf(stderr, message, args);
  fprintf(stderr, "\n");

}
