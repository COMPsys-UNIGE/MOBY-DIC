#include "controller.h"

void scaleP(fxd_in p_in[nP], fxd p_reg[nP])
{
	#pragma HLS ARRAY_PARTITION variable=p_in dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=p_reg dim=1 complete

	fxd_conv p_in_tmp[nP];
	#pragma HLS ARRAY_PARTITION variable=p_in_tmp dim=1 complete

  for (int i = 0; i < nP; i++)
  {
    #pragma HLS PIPELINE
    p_in_tmp[i] = p_in[i];
    p_in_tmp[i] = p_in_tmp[i]*sim_p_scale_gain[i] + sim_p_scale_bias[i];
    p_reg[i] = p_in_tmp[i];
  }
}
