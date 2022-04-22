#include "controller.h"

void scaleD(fxd_in d_in[nD], fxd d_reg[nD])
{
	#pragma HLS ARRAY_PARTITION variable=d_in dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=d_reg dim=1 complete

	fxd_conv d_in_tmp[nD];
	#pragma HLS ARRAY_PARTITION variable=d_in_tmp dim=1 complete

  for (int i = 0; i < nD; i++)
  {
    #pragma HLS PIPELINE
    d_in_tmp[i] = d_in[i];
    d_in_tmp[i] = d_in_tmp[i]*sim_d_scale_gain[i] + sim_d_scale_bias[i];
    d_reg[i] = d_in_tmp[i];
  }
}
