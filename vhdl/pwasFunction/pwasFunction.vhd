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
use work.pwasFunctionPackage.all;

-- Entity declaration

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

-- Architecture declaration

architecture Behavioral of pwasFunction is

	-- Components declaration
	
	component pwasFunctionCircuit is
		Port ( clk : in std_logic;
		     reset : in std_logic;
		     start : in std_logic;
		         x : in x_matrix_pwas;
		     mul_y : in mul_out_matrix_pwas;
		    mul_x1 : out mul_in_matrix_pwas;
		    mul_x2 : out mul_in_matrix_pwas;
		         y : out y_matrix_pwas;
		      done : out std_logic);
	end component;
	
	component mulBank is
		Port (  x1 : in mul_in_matrix_pwas;
		        x2 : in mul_in_matrix_pwas;
		         y : out mul_out_matrix_pwas);
	end component;
	
	-- Signals declaration
	
	signal signal_y: y_matrix_pwas;
	signal signal_x: x_matrix_pwas;
	signal signal_mul_x1: mul_in_matrix_pwas;
	signal signal_mul_x2: mul_in_matrix_pwas;
	signal signal_mul_y: mul_out_matrix_pwas;

begin

	-- Port maps
	
	inst_pwasFunctionCircuit: pwasFunctionCircuit
		port map ( clk => clk,
		         reset => reset,
		         start => start,
		             x => signal_x,
		         mul_y => signal_mul_y,
		        mul_x1 => signal_mul_x1,
		        mul_x2 => signal_mul_x2,
		             y => signal_y,
		          done => done);
					 
	inst_mulBank : mulBank 
	port map (  x1 => signal_mul_x1,
		         x2 => signal_mul_x2,
		          y => signal_mul_y);
			
	-- Assignments

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

end Behavioral;

