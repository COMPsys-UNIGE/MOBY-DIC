----------------------------------------------------------------------------------
-- MOBY-DIC Project
-- www.mobydic-project.eu
--
-- Copyright:
-- (C) 2011 Tomaso Poggi, Spain, tpoggi@essbilbao.org
-- (C) 2011 Alberto Oliveri, University of Genoa, Italy, alberto.oliveri@unige.it
--
-- Legal note:
-- This program is free software; you can redistribute it and/or
-- modify it under the terms of the GNU General Public
-- License as published by the Free Software Foundation; either
-- version 2.1 of the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
-- General Public License for more details.
-- 
-- You should have received a copy of the GNU General Public
-- License along with this library; if not, write to the 
-- Free Software Foundation, Inc., 
-- 59 Temple Place, Suite 330, 
-- Boston, MA  02111-1307  USA
----------------------------------------------------------------------------------

library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

package pwasFunctionPackage is

	-- Declare constants
	
--- BEGIN MATLABGEN ---

	-- Number of domain dimensions
	constant N_DIM_PWAS : integer := 4;
	-- Number of codomain dimensions
	constant N_FUN_PWAS : integer := 1;
	-- Number of bits to represent the inputs
	constant N_BIT_PWAS : integer := 12;
	-- Number of bits to represent the outputs
	constant N_BIT_OUT_PWAS : integer := 12;
	-- Number of bits to represent the coefficients
	constant N_BIT_COEFF_PWAS : integer := 12;
	-- Maximum number of bits used to represent the integer part of the scaled inputs
	constant N_INT_MAX_PWAS : integer := 4;
	-- Number of bits to represent memory address
	constant N_BIT_ADDR_PWAS : integer := 13;
	-- Number of bits for the multipliers
	constant N_BIT_MUL_PWAS : integer := 13;
	
	-- Definition of new types
	
	type x_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(N_BIT_PWAS downto 0);
	type int_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of unsigned(N_INT_MAX_PWAS-1 downto 0);
	type dec_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	type y_matrix_pwas is array(N_FUN_PWAS-1 downto 0) of signed(N_BIT_OUT_PWAS-1 downto 0);
	type mul_in_matrix_pwas is array(N_DIM_PWAS downto 0) of signed(N_BIT_MUL_PWAS-1 downto 0);
	type mul_out_matrix_pwas is array(N_DIM_PWAS downto 0) of signed(2*N_BIT_MUL_PWAS-1 downto 0);
	type mu_matrix_pwas is array(N_DIM_PWAS downto 0) of unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	type coeff_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);
	type coeff_max_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(2*N_BIT_MUL_PWAS-1 downto 0);
	type d_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of std_logic_vector(N_DIM_PWAS-1 downto 0);
	type addr_matrix_pwas is array(N_DIM_PWAS downto 0) of unsigned(N_BIT_ADDR_PWAS-1 downto 0);
	type mem_element_pwas is array(N_DIM_PWAS downto 0) of signed(N_BIT_COEFF_PWAS downto 0);
	type mem_matrix_pwas is array(N_FUN_PWAS-1 downto 0) of mem_element_pwas;
	type integer_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of integer;
	type integer_matrix_out_pwas is array(N_FUN_PWAS-1 downto 0) of integer;

	-- Declare other constants

	-- Number of bits to represent the integer part of the inputs after scaling
	constant N_BIT_INT_PWAS: integer_matrix_pwas := (4,1,4,4);
	constant NP_PWAS: coeff_max_matrix_pwas := ("00010001111111111111111111","00001111111111111111111111","00010001111111111111111111","00010001111111111111111111");
	constant ZERO_PWAS : signed(2*N_BIT_MUL_PWAS-1 downto 0) := (others => '0');
	-- Point position of scaled terms A(x-b)
	constant POINT_POSITION_PWAS : integer_matrix_pwas := (19,22,19,19);
	constant scaleInput_A_PWAS : coeff_matrix_pwas := ("010010000000","010000000000","010010000000","010010000000");
	constant scaleInput_b_PWAS : x_matrix_pwas := ("0000000000000","0000000000000","0000000000000","0000000000000");

	constant START_BIT_PWAS : integer := (11);

--- END MATLABGEN ---

end pwasFunctionPackage;
