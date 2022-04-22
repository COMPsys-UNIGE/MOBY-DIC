#include "controller.h"

// Projects z inside [z_min, z_max]
void project(fxd z[nDim_CTRL]) {

	fxd s1[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=s1 dim=1 complete
	fxd s2[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=s2 dim=1 complete

	for (int i = 0; i < nDim_CTRL; i++)
	{
		#pragma HLS UNROLL
		s1[i] = z_min[i] - z[i];
		s2[i] = z[i] - z_max[i];
	}

	for(int i = 0; i < nDim_CTRL; i++)
	{
		#pragma HLS UNROLL
		if (s1[i] > 0)
			z[i] = z_min[i];
		else if (s2[i] > 0)
			z[i] = z_max[i];
	}
}
