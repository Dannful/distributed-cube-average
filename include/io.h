#pragma once

#include <mpi.h>
#include <stddef.h>

// Point-to-point transfer of float arrays whose element count may exceed
// INT_MAX (the MPI count limit): split into chunks sent in order.
void dc_send_floats(const float *buffer, size_t count, int destination,
                    int tag, MPI_Comm comm);
void dc_recv_floats(float *buffer, size_t count, int source, int tag,
                    MPI_Comm comm);
