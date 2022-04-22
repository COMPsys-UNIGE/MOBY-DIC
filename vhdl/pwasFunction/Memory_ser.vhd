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

entity Memory is
	Port ( clk : in std_logic;
	     start : in std_logic;
	     reset : in std_logic;
	      addr : in addr_matrix_pwas;
	       res : out mem_matrix_pwas;
	      done : out std_logic);
end Memory;

-- Architecture declaration

architecture Behavioral of Memory is
		
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---
    
    proc_state : process(clk,reset)
	begin
		if reset = '0' then
			state <= 0;
		elsif rising_edge(clk) then
			if state = 0 then
				if start = '1' then
					state <= state+1;
				else
					state <= state;
				end if;
			elsif state = N_DIM_PWAS+1 then
				state <= 0;
			else
				state <= state+1;
			end if;
		end if;
	end process;	


	done <= '1' when state = N_DIM_PWAS+1 else
	        '0';
	

end Behavioral;
