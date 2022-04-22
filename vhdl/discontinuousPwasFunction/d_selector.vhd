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

entity d_selector is
	Port ( x_dec_sorted : in dec_matrix_pwas;
	              x_dec : in dec_matrix_pwas;
	                  d : out d_matrix_pwas);
end d_selector;

architecture Behavioral of d_selector is

begin

	process (x_dec, x_dec_sorted)
	begin
		for i in 0 to N_DIM_PWAS-1 loop
			for j in 0 to N_DIM_PWAS-1 loop
				if x_dec(j) < x_dec_sorted(i) then
					d(i)(j) <= '0';
				else 
					d(i)(j) <= '1';
				end if;
			end loop;
		end loop;
	end process;
	
end Behavioral;

