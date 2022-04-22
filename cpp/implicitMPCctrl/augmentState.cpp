#include "controller.h"

void augmentState(fxd x_reg[nX], fxd ref_reg[nX], fxd x_aug[nX_CTRL], fxd ref_aug[nX_CTRL], fxd u_old[nU])
{
  for (int i = 0; i < nX; i++)
  {
      #pragma HLS UNROLL
      x_aug[i] = x_reg[i];
      ref_aug[i] = ref_reg[i];
  }
  for (int i = nX; i < nX_CTRL; i++)
  {
      #pragma HLS UNROLL
      x_aug[i] = u_old[i-nX];
      ref_aug[i] = 0;
  }
}
