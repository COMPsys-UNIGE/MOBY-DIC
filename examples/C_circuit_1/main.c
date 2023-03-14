/* Auto generated code for implementation of embedded system.
   Final user must implement methods configADC(), configDAC()
   readInputADC(y), setOutputDAC(u) and configTimerInterrupt(To)
   as described above:
   configADC() - in this function ADC must be set properly
   configDAC() - in this function DAC must be set properly
   configTimerInterrupt(float interruptTic) - in this function 
         timer and input must be set properly in order to
         activate interrupt every To seconds
   readInputADC(float *y) - in this function the progrma must
         read value of all system output
   setOutputDAC() - in this function the progrma must write
         to DAC controls value
   It is also necessary to call function interruptTimerRoutine()
   into the timer interrupt routine.                            */
#include <stdio.h>
#include "controller.h"
void scaleInput(int xADC[nDim_CTRL], float xscale[nDim_CTRL]);
void scaleOutput(int uADC[nY_CTRL], float uscale[nY_CTRL]);


char controlGo = 0;


void interruptTimerRoutine()
{
controlGo = 1;
}



int main()
{
int xADC[nDim_CTRL]
float x[nDim_CTRL];//input vector
float u[nY_CTRL];//vector containing control value 

int uADC[nY_CTRL]
configADC();//config ADC 
configDAC();//config DAC  

configTimerInterrupt(3.333333e-05);//config timer interrupt to be activate every observer sampling time  

while(1)
{
while(!controlGo);//wait for interrupt activation
controlGo = 0;
readInputADC(xADC);
scaleInput(xADC,x);
/* x = [x1,x2,...,xn,p1,..,pq,d1,..,dz] */

control(x,u);
scaleOutput(uADC,u);
setOutputDAC(uADC);
}
}

void scaleInput(int xADC[nDim_CTRL], float xscale[nDim_CTRL])
{
int i;
for(i=0;i<nDim_CTRL;i++)
xscale[i] = 4.8840048840e+00*xADC[i]+2.4420024420e+00;
xscale[i] = 4.4550605641e+00*xADC[i]+8.8049102528e+02;
xscale[i] = 9.7680097680e-04*xADC[i]+2.0004884005e+00;
xscale[i] = 1.4652014652e-03*xADC[i]+3.0007326007e+00;
}
void scaleOutput(int uADC[nY_CTRL], float uscale[nY_CTRL])
{
int i;
for(i=0;i<nY_CTRL;i++)
uADC[i] = 5.8167675439e+02*uscale[i]+-1.6907465856e+03;
}
