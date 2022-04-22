#include "controller.h"

void admmInit(float x[nX_CTRL], float v1[nDim_CTRL], float v2[nCon + N*nX_CTRL])
{

  float x_reg[nX_CTRL];

  float v2_reg[nCon + N*nX_CTRL];

  float b[N*nX_CTRL];

	float v1_tmp_1[nDim_CTRL];
	float v1_tmp_2[nDim_CTRL];

	float v1_reg_1[nDim_CTRL];
	float v1_reg_2[nDim_CTRL];

/*p
  float p_reg[nP];
  float Ep[nX_CTRL];
p*/

/*d
  float d_reg[nD];
  float Fd[nX_CTRL];
d*/

/*ref
	float ref_reg[nDim_CTRL];

  float q[nDim_CTRL];
  float q_tmp = 0;
  float p_tmp = 0;
ref*/

  // counter
	int k = 1;

  for (int i = 0; i < nX_CTRL; i++)
  {
    x_reg[i] = x[i];

/*ref
    ref_reg[i] = ref[i];
ref*/
  }

  for (int i = 0; i < nCon + N*nX_CTRL; i++)
  {
    v2_reg[i] = 0;
  }

  for (int i = 0; i < N*nX_CTRL; i++)
  {
    b[i] = 0;
  }

/*p
  for (int i = 0; i < nP; i++)
  {
    p_reg[i] = p[i];
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    Ep[i] = 0;
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    for (int j = 0; j < nP; j++)
    {
      Ep[i] += E[i][j]*p_reg[j];
    }
  }
p*/

/*d
  for (int i = 0; i < nD; i++)
  {
    d_reg[i] = d[i];
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    Fp[i] = 0;
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    for (int j = 0; j < nD; j++)
    {
      Fd[i] += F[i][j]*d_reg[j];
    }
  }
d*/

/*ref
  for (int i = 0; i < nDim_CTRL; i++)
  {
    q[i] = 0;
  }
ref*/

  for (int i = 0; i < nX_CTRL; i++)
  {
    for (int j = 0; j < nX_CTRL; j++)
    {
      b[i] += A_pwr[i][j]*x_reg[j];
    }
  }

  for (int i = 0; i < N*nX_CTRL; i+=nX_CTRL)
  {
    for (int j = 0; j < nX_CTRL; j++)
    {

      b[i+j] += G[j];

      /*p
      b[i+j] += Ep[j];
      p*/

      /*d
      b[i+j] += Fd[j];
      d*/
    }
  }

  for (int i = 0; i < nCon; i++)
  {
    v2_reg[i] = v3[i];
  }
  for (int i = nCon; i < nCon + N*nX_CTRL; i++)
  {
    v2_reg[i] = b[i-nCon];
  }

  for (int i = 0; i < nDim_CTRL; i++)
  {
    v1_tmp_1[i] = 0;
    v1_tmp_2[i] = 0;
  }

  for (int j = 0; j < nDim_CTRL; j++)
  {

    // y update
    for (int i = 0; i < nCon + N*nX_CTRL; i++)
    {
      v1_tmp_1[j] -= M2[j][i] * v2_reg[i];
    }
  }

/*ref
  for (int i = 0; i < nX_CTRL; i++)
  {
    q_tmp = 0;
    p_tmp = 0;
    for (int j = 0; j < nX_CTRL; j++)
    {
      q_tmp -= Q[i][j] * ref_reg[j];
    }
    for (int j = 0; j < nX_CTRL; j++)
    {
      p_tmp -= P[i][j] * ref_reg[j];
    }
    for (int l = Nu*nU; l < nDim_CTRL-nX_CTRL; l+=nX_CTRL)
    {
      q[i + l] = 2*q_tmp;
    }
    q[i + nDim_CTRL-nX_CTRL] = 2*p_tmp;
  }
ref*/

  for (int j = 0; j < nDim_CTRL; j++)
  {
    // y update
    for (int i = 0; i < nDim_CTRL; i++)
    {
      v1_tmp_2[j] += M1[j][i] * q[i];
    }
  }

  for (int i = 0; i < nDim_CTRL; i++)
  {
    v1_reg_1[i] = v1_tmp_1[i];
    v1_reg_2[i] = v1_tmp_2[i];
  }

  for (int i = 0; i < nDim_CTRL; i++)
  {
    v1[i] = v1_reg_1[i] + v1_reg_2[i];
  }

  for (int i = 0; i < nCon + N*nX_CTRL; i++)
  {
    v2[i] = v2_reg[i];
  }
}
