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

	component bitonic8 is
		Port ( clk : in std_logic;
	        reset : in std_logic;
	        start : in std_logic;
	           x1 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x2 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x3 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x4 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x5 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x6 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x7 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           x8 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           up : in std_logic;
	           y1 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y2 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y3 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y4 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y5 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y6 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y7 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	           y8 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	         done : out std_logic);
	end component;
	
	signal signal_x1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x3 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x4 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x5 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x6 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x7 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_x8 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y3 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y4 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y5 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y6 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y7 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal signal_y8 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	
begin

	inst_bitonic_8: bitonic8
		port map ( clk => clk,
	            reset => reset,
	            start => start,
	               x1 => signal_x1,
	               x2 => signal_x2,
	               x3 => signal_x3,
	               x4 => signal_x4,
	               x5 => signal_x5,
	               x6 => signal_x6,
	               x7 => signal_x7,
	               x8 => signal_x8,
	               y1 => signal_y1,
	               y2 => signal_y2,
	               y3 => signal_y3,
	               y4 => signal_y4,
	               y5 => signal_y5,
	               y6 => signal_y6,
	               y7 => signal_y7,
	               y8 => signal_y8,
	               up => '1',
	             done => done);

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

end Behavioral;

