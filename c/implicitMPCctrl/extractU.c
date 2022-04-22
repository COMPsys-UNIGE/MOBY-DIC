#include "controller.h"

void extractU(float z[nDim_CTRL], float u_reg[nU], float u_old[nU])
{
  // delta output
  float delta_u[nU];

  for (int i = 0; i < nU; i++)
  {
    if (z[i] < (z_min[Nu*nU+nX+i] - u_old[i]))
        delta_u[i] = z_min[Nu*nU+nX+i] - u_old[i];
    else if (z[i] > (z_max[Nu*nU+nX+i] - u_old[i]))
        delta_u[i] = z_max[Nu*nU+nX+i] - u_old[i];
    else
        delta_u[i] = z[i];
    u_reg[i] = u_old[i] + delta_u[i];
    u_old[i] = u_reg[i];
  }
}
