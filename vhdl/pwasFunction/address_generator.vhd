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

entity address_generator is
	Port (   clk : in std_logic;
	       reset : in std_logic;
	       start : in std_logic;
	       x_int : in int_matrix_pwas;
	       x_dec : in dec_matrix_pwas;
	x_dec_sorted : in dec_matrix_pwas;
	     address : out addr_matrix_pwas;
	        done : out std_logic);
end address_generator;

architecture Behavioral of address_generator is

	-- Components declaration

	component d_selector
		Port ( x_dec_sorted : in  dec_matrix_pwas;
	                 x_dec : in  dec_matrix_pwas;
	                     d : out  d_matrix_pwas);
	end component;

	-- Signals declaration

	signal signal_d : d_matrix_pwas;
	 
begin

	-- Port maps

	inst_d_selector: d_selector 
		port map ( x_dec_sorted => x_dec_sorted,
		                  x_dec => x_dec,
		                      d => signal_d);

	proc_address: process(clk,reset)
	
		-- Stores one bit of signal_d
		variable d_bit: unsigned(0 downto 0);
		variable index: integer;
		
	begin
		
		if reset = '0' then
			for j in 0 to N_DIM_PWAS loop
				address(j) <= (others => '0');
			end loop;
		elsif rising_edge(clk) and start = '1' then
			-- The first address is obtained by concatenating the integer parts
			-- of the input point
			index := -1;
			for j in 0 to N_DIM_PWAS-1 loop
				index := index + N_BIT_INT_PWAS(j);
				address(0)(index downto index-N_BIT_INT_PWAS(j)+1) <= x_int(j)(N_BIT_INT_PWAS(j)-1 downto 0) ;
			end loop;
			
			-- The other addresses are obtained by concatenating the integer parts
			-- of the input points, plus the corresponding value of d
			for i in 1 to N_DIM_PWAS loop           
				index := -1;
				for j in 0 to N_DIM_PWAS-1 loop
					index := index + N_BIT_INT_PWAS(j);
					d_bit := ""&signal_d(i-1)(j);
					address(i)(index downto index-N_BIT_INT_PWAS(j)+1) <= x_int(j)(N_BIT_INT_PWAS(j)-1 downto 0)+d_bit;
				end loop;
			end loop;
		end if;
	end process;
	
	done <= '0' when reset = '0' else
	        start when rising_edge(clk);
	
end Behavioral;

