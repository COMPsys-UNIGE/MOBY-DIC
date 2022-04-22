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
use work.discontinuousPwasFunctionPackage.all;

entity Memory is
	Port ( clk : in std_logic;
	     start : in std_logic;
	     reset : in std_logic;
	      addr : in addr_matrix_pwas;
       dynamic : in integer range 0 to N_DYN_PWAS-1;
	       res : out mem_matrix_pwas;
	      done : out std_logic);
end Memory;

-- Architecture declaration

architecture Behavioral of Memory is
		
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---
	
	done <= '0' when reset = '0' else
	        start when rising_edge(clk);
	

end Behavioral;
