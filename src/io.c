#include "../include/io.h"
#include "mpi.h"

void MPI_Safe_Recv(void *buffer, size_t count, MPI_Datatype data_type, unsigned int source, int tag, MPI_Comm channel) {
  MPI_Status status;
  int received = 0;
  while(received < count) {
    MPI_Recv(buffer + received, count, data_type, source, tag, channel, &status);
    MPI_Get_count(&status, data_type, &received);
  }
}
