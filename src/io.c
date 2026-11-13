#include "io.h"

#define DC_MAX_CHUNK ((size_t)1 << 30)

void dc_send_floats(const float *buffer, size_t count, int destination,
                    int tag, MPI_Comm comm) {
  for (size_t offset = 0; offset < count; offset += DC_MAX_CHUNK) {
    size_t chunk = count - offset < DC_MAX_CHUNK ? count - offset : DC_MAX_CHUNK;
    MPI_Send(buffer + offset, (int)chunk, MPI_FLOAT, destination, tag, comm);
  }
}

void dc_recv_floats(float *buffer, size_t count, int source, int tag,
                    MPI_Comm comm) {
  for (size_t offset = 0; offset < count; offset += DC_MAX_CHUNK) {
    size_t chunk = count - offset < DC_MAX_CHUNK ? count - offset : DC_MAX_CHUNK;
    MPI_Recv(buffer + offset, (int)chunk, MPI_FLOAT, source, tag, comm,
             MPI_STATUS_IGNORE);
  }
}
