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
use IEEE.std_logic_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use work.pwasFunctionPackage.all;

entity Sorter is
	Port ( clk : in std_logic;
	     start : in std_logic;
	     reset : in std_logic;
	         x : in dec_matrix_pwas;
	         y : out dec_matrix_pwas;
	      done : out std_logic);
end Sorter;

architecture Behavioral of Sorter is

	component bitonic2 is
	Port ( x1 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       x2 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       y1 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       y2 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       up : in std_logic;
	    start : in std_logic;
	    reset : in std_logic;
	      clk : in std_logic;
	     done : out std_logic);
	end component;
	
	signal signal_x1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	
begin

	inst_bitonic_2: bitonic2
	port map ( clk => clk,
	         reset => reset,
	         start => start,
	            x1 => signal_x1,
	            x2 => signal_x2,
	            y1 => signal_y1,
	            y2 => signal_y2,
	            up => '1',
	          done => done);
		
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

end Behavioral;

