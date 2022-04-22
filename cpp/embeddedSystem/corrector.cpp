#include "embeddedSystem.h"

void corrector(fxd y[nY], fxd p[nP], fxd u_old[nU], fxd currentState[nX_OBS])
{
    fxd obsInput[nU+nP+nY];

    fxd currentState_reg[nX_OBS];

    fxd_tmp_p predState[nX_OBS];

    for (int i = 0; i < nX_OBS; i++)
    {
      currentState_reg[i] = currentState[i];
    }

    for (int i = 0; i < nU; i++)
    {
      obsInput[i] = u_old[i];
    }

    /*p
    for (int i = 0; i < nP; i++)
    {
      obsInput[nU+i] = p[i];
    }
    p*/

    for (int i = 0; i < nY; i++)
    {
      obsInput[nU+nP+i] = y[i];
    }

    // Update current state
    for(int i = 0; i < nX_OBS; i++)
    {
      #pragma HLS UNROLL //small off
      predState[i] = 0;
      for(int j = 0; j < nX_OBS; j++)
      {
        #pragma HLS UNROLL //small
        predState[i] += Fupdate[i][j]*currentState_reg[j];
      }
      for(int j = 0; j < nU + nP + nY; j++)
      {
        #pragma HLS UNROLL //small
        predState[i] += Fupdate[i][j+nX_OBS]*obsInput[j];
      }
      predState[i] += Gupdate[i];
    }

    for (int i = 0; i < nX_OBS; i++)
    {
      currentState[i] = predState[i];
    }
}
