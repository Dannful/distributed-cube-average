#include "../include/io.h"
#include <mpi.h>

void MPI_Safe_Recv(void *buffer, size_t count, MPI_Datatype data_type, unsigned int source, int tag, MPI_Comm channel) {
  MPI_Status status;
  int received_total = 0;
  int data_type_size;
  MPI_Type_size(data_type, &data_type_size);
  while(received_total < count) {
    int received;
    MPI_Recv((char*) buffer + received_total * data_type_size, count - received_total, data_type, source, tag, channel, &status);
    MPI_Get_count(&status, data_type, &received);
    received_total += received;
  }
}
