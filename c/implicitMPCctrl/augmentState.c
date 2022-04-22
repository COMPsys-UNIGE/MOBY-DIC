#include "controller.h"

void augmentState(float x_reg[nX], float ref_reg[nX], float x_aug[nX_CTRL], float ref_aug[nX_CTRL], float u_old[nU])
{
  for (int i = 0; i < nX; i++)
  {
      x_aug[i] = x_reg[i];
      ref_aug[i] = ref_reg[i];
  }
  for (int i = nX; i < nX_CTRL; i++)
  {
      x_aug[i] = u_old[i-nX];
      ref_aug[i] = 0;
  }
}
