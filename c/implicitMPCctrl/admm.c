#include "controller.h"

// Alternating Direction Method of Multipliers (ADMM)
// x is the current system state, q is a reference dependent vector;
// z and lambda are initial values of the optimization variables and
// Lagrange multipliers, respectively
//void admm(float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL], float z[nDim_CTRL], float lambda[nDim_CTRL])
void admm(float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL], float z[nDim_CTRL])
{

	// internal copy of y
	static float y[nCon + N*nX_CTRL];

	// internal copy of lambda
	static float lambda[nCon + N*nX_CTRL];

	// internal copy of b
	float v1_reg[nDim_CTRL];

	float v2_reg[nCon + N*nX_CTRL];

	float z_reg[nDim_CTRL];

	float z_tmp[nDim_CTRL];
	float y_tmp[nDim_CTRL];
	float lambda_tmp[nCon + N*nX_CTRL];

	float lambda_new[nCon + N*nX_CTRL];

	static bool init = true;

	if (init == true)
	{
			for (int i = 0; i < nCon + N*nX_CTRL; i++)
			{
				y[i] = 0;
				lambda[i] = 0;
			}

			init = false;
	}

	for (int i = 0; i < nCon + N*nX_CTRL; i++)
	{
		v2_reg[i] = v2[i];
	}

	for (int i = 0; i < nDim_CTRL; i++)
	{
		v1_reg[i] = v1[i];
	}

	// ADMM iterations
	for (int k = 0; k < max_iter; k++)
	{
		for (int i = 0; i < nDim_CTRL; i++)
		{
			z_tmp[i] = 0;
		}

		for (int j = 0; j < nDim_CTRL; j++)
		{
			// y update
			for (int i = 0; i < nCon + N*nX_CTRL; i++)
			{
				z_tmp[j] += M2[j][i] * (y[i] + lambda[i]);
			}
		}

		for (int i = 0; i < nDim_CTRL; i++)
		{
			z_reg[i] = z_tmp[i];
			z_reg[i] += v1_reg[i];
		}

		// z update (projection)
		for (int i = 0; i < nDim_CTRL; i++)
		{
			y_tmp[i] = 0;
		}

		for (int j = 0; j < nCon; j++)
		{
			// y update
			for (int i = 0; i < nDim_CTRL; i++)
			{
				y_tmp[j] += M3[j][i] * z_reg[i];
			}
		}

		for (int i = 0; i < nCon; i++)
		{
			y[i] = y_tmp[i];
			y[i] += v3[i] - lambda[i];

			if (y[i] < 0)
			{
				y[i] = 0;
			}
		}

		// lambda update
		for (int i = 0; i < nCon + N*nX_CTRL; i++)
		{
			lambda_tmp[i] = 0;
		}

		for (int j = 0; j < nCon + N*nX_CTRL; j++)
		{
			// y update
			for (int i = 0; i < nDim_CTRL; i++)
			{
				lambda_tmp[j] += M4[j][i] * z_reg[i];
			}
		}

		for (int i = 0; i < nCon + N*nX_CTRL; i++)
		{
			lambda_new[i] = lambda_tmp[i];
		}

		for (int i = 0; i < nCon + N*nX_CTRL; i++)
		{
			lambda[i] += lambda_new[i] + y[i] - v2_reg[i];
		}

	}

	for (int i = 0; i < nDim_CTRL; i++)
	{
		z[i] = z_reg[i];
	}
}
