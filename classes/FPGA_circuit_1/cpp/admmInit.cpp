#include "controller.h"

//void admmInit(fxd x[nX_CTRL], fxd z[nDim_CTRL], fxd lambda[nDim_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + (N+1)*nX_CTRL])
void admmInit(fxd x[nX_CTRL], fxd ref[nX_CTRL], fxd v1[nDim_CTRL], fxd v2[nCon + (N+1)*nX_CTRL])
{
  #pragma HLS ARRAY_PARTITION variable=x dim=1 complete
  #pragma HLS ARRAY_PARTITION variable=v1 dim=1 complete
  #pragma HLS ARRAY_PARTITION variable=v2 dim=1 complete

  /*p
  #pragma HLS ARRAY_PARTITION variable=p dim=1 complete
  p*/

  /*d
  #pragma HLS ARRAY_PARTITION variable=d dim=1 complete
  d*/

  
  #pragma HLS ARRAY_PARTITION variable=ref dim=1 complete
  

  //static bool firstCall = true;

  fxd x_reg[nX_CTRL];
  #pragma HLS ARRAY_PARTITION variable=x_reg dim=1 complete

  fxd v2_reg[nCon + (N+1)*nX_CTRL];
  #pragma HLS ARRAY_PARTITION variable=v2_reg dim=1 complete

  fxd b[(N+1)*nX_CTRL];
  #pragma HLS ARRAY_PARTITION variable=b dim=1 complete

	fxd_tmp_m v1_tmp_1[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=v1_tmp_1 dim=1 complete
	fxd_tmp_m v1_tmp_2[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=v1_tmp_2 dim=1 complete
	fxd_tmp_m v1_tmp[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=v1_tmp dim=1 complete

	// fxd v1_reg_1[nDim_CTRL];
	// #pragma HLS ARRAY_PARTITION variable=v1_reg_1 dim=1 complete
	// fxd v1_reg_2[nDim_CTRL];
	// #pragma HLS ARRAY_PARTITION variable=v1_reg_2 dim=1 complete

/*p
  fxd p_reg[nP];
  #pragma HLS ARRAY_PARTITION variable=p_reg dim=1 complete
  fxd Ep[nX_CTRL];
  #pragma HLS ARRAY_PARTITION variable=Ep dim=1 complete
p*/

/*d
  fxd d_reg[nD];
  #pragma HLS ARRAY_PARTITION variable=d_reg dim=1 complete
  fxd Fd[nX_CTRL];
  #pragma HLS ARRAY_PARTITION variable=Fd dim=1 complete
d*/


	fxd ref_reg[nDim_CTRL];
	#pragma HLS ARRAY_PARTITION variable=ref_reg dim=1 complete

  fxd q[nDim_CTRL];
  #pragma HLS ARRAY_PARTITION variable=q dim=1 complete
  fxd q_tmp = 0;
  fxd p_tmp = 0;


  // counter
	int k = 1;

  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    x_reg[i] = x[i];


    ref_reg[i] = ref[i];

  }

  for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    v2_reg[i] = 0;
  }

  for (int i = 0; i < (N+1)*nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    b[i] = 0;
  }

/*p
  for (int i = 0; i < nP; i++)
  {
    #pragma HLS UNROLL
    p_reg[i] = p[i];
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    Ep[i] = 0;
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS PIPELINE off
    for (int j = 0; j < nP; j++)
    {
      #pragma HLS PIPELINE
      Ep[i] += E[i][j]*p_reg[j];
    }
  }
p*/

/*d
  for (int i = 0; i < nD; i++)
  {
    #pragma HLS UNROLL
    d_reg[i] = d[i];
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    Fp[i] = 0;
  }

  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS PIPELINE off
    for (int j = 0; j < nD; j++)
    {
      #pragma HLS PIPELINE
      Fd[i] += F[i][j]*d_reg[j];
    }
  }
d*/


  for (int i = 0; i < nDim_CTRL; i++)
  {
    #pragma HLS UNROLL
    q[i] = 0;
  }


  // for (int i = 0; i < nX_CTRL; i++)
  // {
  //   #pragma HLS PIPELINE off
  //   for (int j = 0; j < nX_CTRL; j++)
  //   {
  //     #pragma HLS PIPELINE
  //     b[i] += A_pwr[i][j]*x_reg[j];
  //   }
  // }

  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS UNROLL
      b[i] = x_reg[i];
  }

  for (int i = nX_CTRL; i < (N+1)*nX_CTRL; i+=nX_CTRL)
  {
    //#pragma HLS PIPELINE II=6
    //#pragma HLS PIPELINE off
    for (int j = 0; j < nX_CTRL; j++)
    {
      #pragma HLS PIPELINE

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
    #pragma HLS UNROLL
    v2_reg[i] = v3[i];
  }
  for (int i = nCon; i < nCon + (N+1)*nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    v2_reg[i] = b[i-nCon];
  }

  for (int i = 0; i < nDim_CTRL; i++)
  {
    #pragma HLS UNROLL
    v1_tmp_1[i] = 0;
    v1_tmp_2[i] = 0;
  }

  for (int j = 0; j < nDim_CTRL; j++)
  {
    #pragma HLS PIPELINE off

    // y update
    for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
    {
      #pragma HLS PIPELINE
      v1_tmp_1[j] -= M2[j][i] * v2_reg[i];
    }
  }

  //for (int i = 0; i < nDim_CTRL; i++)
  //{
    //#pragma HLS UNROLL
    //v1_reg_1[i] = v1_tmp_1[i];
  //}


  for (int i = 0; i < nX_CTRL; i++)
  {
    #pragma HLS PIPELINE off
    q_tmp = 0;
    p_tmp = 0;
    for (int j = 0; j < nX_CTRL; j++)
    {
      #pragma HLS PIPELINE
      q_tmp -= Q[i][j] * ref_reg[j];
    }
    for (int j = 0; j < nX_CTRL; j++)
    {
      #pragma HLS PIPELINE
      p_tmp -= P[i][j] * ref_reg[j];
    }
    for (int l = Nu*nU; l < nDim_CTRL-nX_CTRL; l+=nX_CTRL)
    {
      #pragma HLS PIPELINE
      q[i + l] = 2*q_tmp;
    }
    q[i + nDim_CTRL-nX_CTRL] = 2*p_tmp;
  }


  //for (int i = 0; i < nDim_CTRL; i++)
  //{
    //#pragma HLS UNROLL
    //v1_tmp[i] = 0;
  //}

  for (int j = 0; j < nDim_CTRL; j++)
  {
    #pragma HLS PIPELINE off

    // y update
    for (int i = 0; i < nDim_CTRL; i++)
    {
      #pragma HLS PIPELINE
      v1_tmp_2[j] += M1[j][i] * q[i];
    }
  }

  for (int i = 0; i < nDim_CTRL; i++)
  {
    #pragma HLS UNROLL
    // v1_reg_1[i] = v1_tmp_1[i];
    // v1_reg_2[i] = v1_tmp_2[i];
    v1_tmp[i] = v1_tmp_1[i] + v1_tmp_2[i];
    v1[i] = v1_tmp[i];
  }

  // for (int i = 0; i < nDim_CTRL; i++)
  // {
  //   #pragma HLS PIPELINE
  //   v1[i] = v1_reg_1[i] + v1_reg_2[i];
  // }

  for (int i = 0; i < nCon + (N+1)*nX_CTRL; i++)
  {
    #pragma HLS UNROLL
    v2[i] = v2_reg[i];
  }
}
