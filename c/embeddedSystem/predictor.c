#include "embeddedSystem.h"

void predictor(float y[nY], float p[nP], float u_old[nU], float currentState[nX_OBS])
{
    float obsInput[nU+nP+nY];

    float currentState_reg[nX_OBS];

    float predState[nX_OBS];

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

    // Predict next state
    for(int i = 0; i < nX_OBS; i++)
    {
      predState[i] = 0;
      for(int j = 0; j < nX_OBS; j++)
      {
        predState[i] += Fpred[i][j]*currentState_reg[j];
      }
      for(int j = 0; j < nU + nP + nY; j++)
      {
        predState[i] += Fpred[i][j+nX_OBS]*obsInput[j];
      }
      predState[i] += Gpred[i];
    }

    for (int i = 0; i < nX_OBS; i++)
    {
      currentState[i] = predState[i];
    }
}
