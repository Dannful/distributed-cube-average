#pragma once

#include <mpi.h>
#include <stddef.h>

void MPI_Safe_Recv(void *buffer, size_t count, MPI_Datatype data_type, unsigned int receiver, int tag, MPI_Comm channel);
