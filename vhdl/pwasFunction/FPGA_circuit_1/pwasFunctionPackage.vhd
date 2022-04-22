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
	constant N_DIM_PWAS : integer := 2;
	-- Number of codomain dimensions
	constant N_FUN_PWAS : integer := 1;
	-- Number of bits to represent the inputs
	constant N_BIT_PWAS : integer := 12;
	-- Number of bits to represent the outputs
	constant N_BIT_OUT_PWAS : integer := 12;
	-- Number of bits to represent the coefficients
	constant N_BIT_COEFF_PWAS : integer := 12;
	-- Maximum number of bits used to represent the integer part of the scaled inputs
	constant N_INT_MAX_PWAS : integer := 3;
	-- Number of bits to represent memory address
	constant N_BIT_ADDR_PWAS : integer := 6;
	-- Number of bits for the multipliers
	constant N_BIT_MUL_PWAS : integer := 13;
	
	-- Definition of new types
	
	type x_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of signed(N_BIT_PWAS-1 downto 0);
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
	type mem_element_pwas is array(N_DIM_PWAS downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);
	type mem_matrix_pwas is array(N_FUN_PWAS-1 downto 0) of mem_element_pwas;
	type integer_matrix_pwas is array(N_DIM_PWAS-1 downto 0) of integer;
	type integer_matrix_out_pwas is array(N_FUN_PWAS-1 downto 0) of integer;
	type p0_matrix_pwas is array(6 downto 0) of signed(N_BIT_PWAS-1 downto 0);
	type p1_matrix_pwas is array(6 downto 0) of signed(N_BIT_PWAS-1 downto 0);
	type scale_matrix0_a_pwas is array(6 downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);
	type scale_matrix0_b_pwas is array(6 downto 0) of signed(N_BIT_PWAS-1 downto 0);
	type scale_matrix1_a_pwas is array(6 downto 0) of signed(N_BIT_COEFF_PWAS-1 downto 0);
	type scale_matrix1_b_pwas is array(6 downto 0) of signed(N_BIT_PWAS-1 downto 0);

	-- Declare other constants

	-- Number of bits to represent the integer part of the inputs after scaling
	constant N_BIT_INT_PWAS: integer_matrix_pwas := (3,3);
	constant NP_PWAS: coeff_max_matrix_pwas := ("00000011111111111111111111","00000011111111111111111111");
	constant ZERO_PWAS : signed(2*N_BIT_MUL_PWAS-1 downto 0) := (others => '0');
	-- Point position of scaled terms A(x-b)
	constant POINT_POSITION_PWAS : integer_matrix_pwas := (20,20);
	constant NUMPART_PWAS: integer_matrix_pwas := (7,7);
	constant scaleInput0_A_PWAS : scale_matrix0_a_pwas := ("011100000000","011100000000","011100000000","011100000000","011100000000","011100000000","011100000000");
	constant scaleInput1_A_PWAS : scale_matrix1_a_pwas := ("011100000000","011100000000","011100000000","011100000000","011100000000","011100000000","011100000000");
	constant scaleInput0_b_PWAS : scale_matrix0_b_pwas := ("010110110110","001101101101","000100100100","111011011011","110010010010","101001001001","100000000000");
	constant scaleInput1_b_PWAS : scale_matrix1_b_pwas := ("010110110110","001101101101","000100100100","111011011011","110010010010","101001001001","100000000000");

	constant START_BIT_PWAS : integer := (11);

	constant P0_PWAS : p0_matrix_pwas := ("011111111111","010110110110","001101101101","000100100100","111011011011","110010010010","101001001001");
	constant P1_PWAS : p1_matrix_pwas := ("011111111111","010110110110","001101101101","000100100100","111011011011","110010010010","101001001001");
--- END MATLABGEN ---

end pwasFunctionPackage;
