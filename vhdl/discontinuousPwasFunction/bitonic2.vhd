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
use work.discontinuousPwasFunctionPackage.all;

entity bitonic2 is
	Port ( x1 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       x2 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       y1 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       y2 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	       up : in std_logic;
	    start : in std_logic;
	    reset : in std_logic;
	      clk : in std_logic;
	     done : out std_logic);
end bitonic2;

architecture Behavioral of bitonic2 is

	component comparator is
		Port ( clk : in std_logic;
		     start : in std_logic;
		     reset : in std_logic;
		        x1 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        x2 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        y1 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        y2 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        up : in std_logic);
	end component;

begin

	comp1: comparator
		port map( clk => clk,
		        start => start,
		        reset => reset,
		           x1 => x1,
		           x2 => x2,
		           y1 => y1,
		           y2 => y2,
		           up => up);

	done <= '0' when reset = '0' else
	        start when rising_edge(clk);

end Behavioral;

