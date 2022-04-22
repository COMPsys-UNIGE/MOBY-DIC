# MOBY-DIC Toolbox

**MOBY-DIC Toolbox** is a MATLAB toolbox for the automatic generation of embedded control systems based on Model Predictive Control (MPC).

The main feature of the toolbox is the automatic generation of C and/or VHDL code, for the circuit implementation of controllers (either exact/approximate explicit or implicit) and observers on microcontrollers and/or Field Programmable Gate Arrays (FPGA), respectively. Also, the design of model-based state observers such as Kalman filters or Kalman predictors is included in the toolbox.

Simulink models for the simulation of the whole closed-loop system (comprising controller and observer) and Xilinx System Generator models for hardware-in-the-loop simulations can be easily generated.

The functionalities of **MOBY-DIC Toolbox** are:

- Open-loop simulation (through MATLAB and Simulink) of Linear Time Invariant (LTI) systems comprising measurable parameters and unmeasurable inputs (disturbances)
- Open-loop simulation (through MATLAB and Simulink) of Piecewise Affine (PWA) systems comprising measurable parameters and unmeasurable inputs (disturbances)
- Design of implicit/explicit MPC controllers (for state/output regulation or tracking) for constrained LTI systems
- Design of explicit switched MPC controllers (for state/output regulation or tracking) for constrained PWA systems
- Design of observers (Kalman filters or Kalman predictors) for the state estimation of LTI systems
- Design of switched Kalman filters or Kalman predictors for the state estimation of PWA systems
- Closed-loop simulation (through MATLAB and Simulink) of LTI or PWA systems, including Kalman filter or Kalman predictor
- Hardware-in-the-loop simulation (through Xilinx System Generator) of LTI or PWA systems, including Kalman filter or Kalman predictor
- Approximation of arbitrary nonlinear functions through PWAS functions
- Approximation of explicit MPC control laws through PWAS functions
- Generation of C code for the implementation on microcontroller of PWAS functions
- Generation of VHDL code for the implementation on FPGA of PWAS functions
- Generation of C code for the implementation on microcontroller of controller (either exact/approximate explicit or implicit) and observer
- Generation of VHDL code for the implementation on FPGA of controller (either exact/approximate explicit or implicit) and observer

## Prerequisites
- MATLAB R2011a or later
- Control System Toolbox
- Optimization Toolbox
- Simulink (needed only to perform simulations)
- Multi Parametric Toolbox 3 (only for explicit MPC, downloadable at https://www.mpt3.org/)
- Xilinx Vitis HLS (only for generating VHDL code of implicit MPC)
- A supported compiler (needed only to recompile mex files, see below)

## Installation procedure
Launch script *installMOBYDIC.m* and follow the instructions.
