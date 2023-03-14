#include "controller.h"

// Alternating Direction Method of Multipliers (ADMM)
// x is the current system state, q is a reference dependent vector;
// z and lambda are initial values of the optimization variables and
// Lagrange multipliers, respectively
//void admm(fxd v1[nDim_CTRL], fxd v2[nCon + (N+1)*nX_CTRL], fxd z[nDim_CTRL], fxd lambda[nDim_CTRL])
void admm(fxd v1[nDim_CTRL], fxd v2[nCon + (N+1)*nX_CTRL], fxd z[nDim_CTRL])
{
	#pragma HLS allocation operation instances=mul limit=59
	#pragma HLS ARRAY_PARTITION variable=v1 dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=v2 dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=z dim=1 complete

	// internal copy of y
	static fxd y[nCon + (N+1)*nX_CTRL];
	#pragma HLS ARRAY_PARTITION variable=y dim=1 complete

	// internal copy of lambda
	static fxd lambda[nCon + (N+1)*nX_CTRL];
	#pragma HLS ARRAY_PARTITION variable=lambda dim=1 complete

	// internal copy of b
	fxd v1_reg[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=v1_reg dim=1 complete

	fxd v2_reg[nCon + (N+1)*nX_CTRL];
	#pragma HLS ARRAY_PARTITION variable=v2_reg dim=1 complete

	fxd z_reg[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=z_reg dim=1 complete

	fxd_tmp_m z_tmp[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=z_tmp dim=1 complete
	fxd_tmp_m y_tmp[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=y_tmp dim=1 complete
	fxd_tmp_m lambda_tmp[nCon + (N+1)*nX_CTRL];
	#pragma HLS ARRAY_PARTITION variable=lambda_tmp dim=1 complete

	fxd lambda_new[nCon + (N+1)*nX_CTRL];
	#pragma HLS ARRAY_PARTITION variable=lambda_new dim=1 complete

	static bool init = true;

	if (init == true)
	{
			for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
			{
				#pragma HLS UNROLL
				y[i] = 0;
				lambda[i] = 0;
			}

			init = false;
	}

	for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
	{
		#pragma HLS UNROLL
		v2_reg[i] = v2[i];
	}

	for (int i = 0; i < nDim_CTRL; i++)
	{
		#pragma HLS UNROLL
		v1_reg[i] = v1[i];
	}

	// ADMM iterations
	for (int k = 0; k < max_iter; k++)
	{
		#pragma HLS PIPELINE off

		for (int i = 0; i < nDim_CTRL; i++)
		{
			#pragma HLS UNROLL
			z_tmp[i] = 0;
		}

		for (int j = 0; j < nDim_CTRL; j++)
		{
			//#pragma HLS allocation operation instances=mul limit=59
			#pragma HLS UNROLL //small off

			// y update
			for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
			{
				#pragma HLS UNROLL //small
				z_tmp[j] += M2[j][i] * (y[i] + lambda[i]);
			}
		}

		for (int i = 0; i < nDim_CTRL; i++)
		{
			#pragma HLS UNROLL //small
			z_reg[i] = z_tmp[i];
			z_reg[i] += v1_reg[i];
		}

		// z update (projection)
		for (int i = 0; i < nDim_CTRL; i++)
		{
			#pragma HLS UNROLL
			y_tmp[i] = 0;
		}

		for (int j = 0; j < nCon; j++)
		{
			//#pragma HLS allocation operation instances=mul limit=59
			#pragma HLS UNROLL //small off

			// y update
			for (int i = 0; i < nDim_CTRL; i++)
			{
				#pragma HLS UNROLL //small
				y_tmp[j] += M3[j][i] * z_reg[i];
			}
		}

		for (int i = 0; i < nCon; i++)
		{
			#pragma HLS UNROLL //small
			y[i] = y_tmp[i];
			y[i] += v3[i] - lambda[i];

			if (y[i] < 0)
			{
				y[i] = 0;
			}
		}

		// lambda update
		for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
		{
			#pragma HLS UNROLL
			lambda_tmp[i] = 0;
		}

		for (int j = 0; j < nCon + (N+1)*nX_CTRL; j++)
		{
			//#pragma HLS allocation operation instances=mul limit=59
			#pragma HLS UNROLL //small off

			// y update
			for (int i = 0; i < nDim_CTRL; i++)
			{
				#pragma HLS UNROLL //small
				lambda_tmp[j] += M4[j][i] * z_reg[i];
			}
		}

		for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
		{
			#pragma HLS UNROLL
			lambda_new[i] = lambda_tmp[i];
		}

		for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
		{
			#pragma HLS UNROLL //small
			lambda[i] += lambda_new[i] + y[i] - v2_reg[i];
		}

	}

	for (int i = 0; i < nDim_CTRL; i++)
	{
		#pragma HLS UNROLL
		z[i] = z_reg[i];
	}
}
