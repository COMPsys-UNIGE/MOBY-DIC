#include "controller.h"

void scaleRef(fxd_in ref_in[nRef], fxd ref_reg[nX])
{
  fxd_conv ref_tmp[nX];
	#pragma HLS ARRAY_PARTITION variable=ref_tmp dim=1 complete

  for (int i = 0; i < nX; i++)
  {
    #pragma HLS PIPELINE
  	ref_reg[i] = 0;
  }

	for (int i = 0; i < nRef; i++)
  {
    #pragma HLS PIPELINE
    ref_tmp[ref_idx[i]] = ref_in[i];
    ref_tmp[ref_idx[i]] = ref_tmp[ref_idx[i]]*sim_xref_scale_gain[i] + sim_xref_scale_bias[i];
  	ref_reg[ref_idx[i]] = ref_tmp[ref_idx[i]];
  }
}
