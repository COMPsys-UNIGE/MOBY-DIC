#include "controller.h"

void extractU(fxd z[nDim_CTRL], fxd u_reg[nU], fxd u_old[nU])
{
  // delta output
  fxd delta_u[nU];
  #pragma HLS ARRAY_PARTITION variable=delta_u dim=1 complete

  for (int i = 0; i < nU; i++)
  {
    #pragma HLS UNROLL
    delta_u[i] = z[i];
    u_reg[i] = u_old[i] + delta_u[i];
    u_old[i] = u_reg[i];
  }
}
