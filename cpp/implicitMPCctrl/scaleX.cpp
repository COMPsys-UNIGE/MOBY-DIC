#include "controller.h"

void scaleX(fxd_in x_in[nX], fxd x_reg[nX])
{
	fxd_conv x_in_tmp[nX];
	#pragma HLS ARRAY_PARTITION variable=x_in_tmp dim=1 complete

  for (int i = 0; i < nX; i++)
  {
    #pragma HLS PIPELINE
    x_in_tmp[i] = x_in[i];
    x_in_tmp[i] = x_in_tmp[i]*sim_x_scale_gain[i] + sim_x_scale_bias[i];
    x_reg[i] = x_in_tmp[i];
  }
}
