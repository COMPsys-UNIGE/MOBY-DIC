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

entity Sum is
	Port ( clk : in  std_logic;
	     reset : in  std_logic;
	     start : in  std_logic;
	         x : in  mul_out_matrix_pwas;
	         y : out y_matrix_pwas;
	      done : out  std_logic);
end Sum;

architecture Behavioral of Sum is

	signal y_tmp: y_matrix_pwas;
	signal state : integer range 0 to N_FUN_PWAS;
	signal sat_state : integer range 0 to N_FUN_PWAS-1;

begin

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
			elsif state = N_FUN_PWAS then
				state <= 0;
			else
				state <= state+1;
			end if;
		end if;
	end process;

	proc_sum: process(clk,reset,start)
	
	variable sum : signed(2*N_BIT_MUL_PWAS+N_DIM_PWAS-1 downto 0);
	variable idx : integer range 0 to N_FUN_PWAS-1;
	
	begin
		if reset = '0' then
			for i in 0 to N_FUN_PWAS-1 loop
				y_tmp(i) <= (others => '0');
			end loop;
			idx := 0;
		elsif rising_edge(clk) and start = '1' then
			sum := (others => '0');
			for i in 0 to N_DIM_PWAS loop
				sum := sum+x(i);
			end loop;
			y_tmp(sat_state) <= sum(START_BIT_PWAS(sat_state)+N_BIT_OUT_PWAS-1 downto START_BIT_PWAS(sat_state));
		end if;
	end process;
	
	sat_state <= state when state < N_FUN_PWAS else
	             N_FUN_PWAS-1;
		
	done <= '1' when state = N_FUN_PWAS else
	        '0';

	proc_y: process(clk,reset,state,y_tmp)
	
	begin
		if reset = '0' then
			for i in 0 to N_FUN_PWAS-1 loop
				y(i) <= (others => '0');
			end loop;
		else
			for i in 0 to N_FUN_PWAS-2 loop
				if rising_edge(clk) and state = N_FUN_PWAS-1 then
					y(i) <= y_tmp(i);
				end if;
			end loop;
			y(N_FUN_PWAS-1) <= y_tmp(N_FUN_PWAS-1);
		end if;
	end process;


end Behavioral;

