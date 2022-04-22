#include "embeddedSystem.h"

void scaleY(fxd_in y_in[nY], fxd y_reg[nY])
{
	fxd_conv y_in_tmp[nY];
	#pragma HLS ARRAY_PARTITION variable=y_in_tmp dim=1 complete

	for (int i = 0; i < nY; i++)
	{
		#pragma HLS PIPELINE
    y_in_tmp[i] = y_in[i];
    y_in_tmp[i] = y_in_tmp[i]*sim_y_scale_gain[i] + sim_y_scale_bias[i];
    y_reg[i] = y_in_tmp[i];
	}
}
