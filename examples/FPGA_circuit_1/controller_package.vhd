--
--	Package File Template
--
--	Purpose: This package defines supplemental types, subtypes, 
--		 constants, and functions 
--
--   To use any of the example code shown below, uncomment the lines and modify as necessary
--

library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

package controller_package is

	-- Declare constants
	
--- BEGIN MATLABGEN ---

	-- Number of domain dimensions
	constant N_DIM_CTRL : integer := 4;
	-- Number of codomain dimensions
	constant N_FUN_CTRL : integer := 1;
	-- Number of bits to represent the inputs
	constant N_BIT_CTRL : integer := 12;
	-- Number of bits to represent the outputs
	constant N_BIT_OUT_CTRL : integer := 12;
	-- Number of bits to represent the coefficients
	constant N_BIT_COEFF_CTRL : integer := 12;
	-- Number of bits to represent memory address
	constant N_BIT_MUL_CTRL : integer := 13;
	-- Definition of new types
	
constant sampling_latency_CTRL : integer := 1666;

	constant nx_CTRL : integer := 2;
	constant npar_CTRL : integer := 1;
	constant nd_CTRL : integer := 0;
	constant nu_CTRL : integer := 1;
	constant nxref_CTRL : integer := 1;
	type xType_CTRL is array(nx_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);
	type uType_CTRL is array(nu_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);
	type pType_CTRL is array(npar_CTRL-1 downto 0) of std_logic_vector(N_BIT_CTRL-1 downto 0);
	type x_matrix_CTRL is array(N_DIM_CTRL-1 downto 0) of signed(N_BIT_CTRL-1 downto 0);
	type u_matrix_CTRL is array(N_FUN_CTRL-1 downto 0) of signed(N_BIT_OUT_CTRL-1 downto 0);
	type mul_in_matrix_CTRL is array(N_DIM_CTRL downto 0) of signed(N_BIT_MUL_CTRL-1 downto 0);
	type mul_out_matrix_CTRL is array(N_DIM_CTRL downto 0) of signed(2*N_BIT_MUL_CTRL-1 downto 0);
	type u_matrix_CTRL_toy is array(1 downto 0) of signed(N_BIT_OUT_CTRL-1 downto 0);
	constant defaultOutput : u_matrix_CTRL_toy := ("000101100101","000101100101");

--- END MATLABGEN ---

end controller_package;
