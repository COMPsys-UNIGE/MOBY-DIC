open_project HLSproject
set_top control
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/controller.h
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/controller.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/scaleX.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/scaleRef.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/scaleU.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/augmentState.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/extractU.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/admmInit.cpp
add_files /home/ale/Code/MOBY-DIC/classes/FPGA_circuit_1/cpp/admm.cpp
open_solution "solution1" -flow_target vivado
set_part {xc7z020clg484-1}
create_clock -period 10000 -name default
config_rtl -reset control -reset_level low
csynth_design
close_project
quit
