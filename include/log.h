#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void dc_log_info(int rank, char message[], ...);
void dc_log_error(int rank, char message[], ...);

#ifdef __cplusplus
}
#endif