#include "embeddedSystem.h"

void scaleDout(fxd d_reg[nD], fxd_out_state d_est[nD])
{
	fxd_conv d_tmp[nD];
	#pragma HLS ARRAY_PARTITION variable=d_tmp dim=1 complete

	for (int i = 0; i < nD; i++)
	{
		#pragma HLS PIPELINE
		d_tmp[i] = d_reg[i];
		d_tmp[i] = (d_tmp[i] - sim_d_scale_bias[i])*sim_d_scale_gain[i];
		d_est[i] = d_tmp[i];
	}
}
