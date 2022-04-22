#include "embeddedSystem.h"

void scaleXout(fxd x_reg[nX], fxd_out_state x_est[nX])
{
		fxd_conv x_tmp[nX];
		#pragma HLS ARRAY_PARTITION variable=x_tmp dim=1 complete

		for (int i = 0; i < nX; i++)
		{
				#pragma HLS PIPELINE
				x_tmp[i] = x_reg[i];
				x_tmp[i] = (x_tmp[i] - sim_x_scale_bias[i])*sim_x_scale_gain[i];
				x_est[i] = x_tmp[i];
		}
}
