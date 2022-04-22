#include "controller.h"

void scaleU(fxd u_reg[nU], fxd_out u_opt[nU])
{
	fxd_conv u_tmp[nU];
	#pragma HLS ARRAY_PARTITION variable=u_tmp dim=1 complete

	for (int i = 0; i < nU; i++)
	{
    #pragma HLS PIPELINE
		u_tmp[i] = u_reg[i];
		u_tmp[i] = (u_tmp[i] - sim_u_scale_bias[i])*sim_u_scale_gain[i];
		u_opt[i] = u_tmp[i];
	}
}
