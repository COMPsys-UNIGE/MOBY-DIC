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

entity mulBankManager is
	Port (     clk : in std_logic;
	         reset : in std_logic;
	         start : in std_logic;
	   sorter_done : in std_logic;
	      mem_done : in std_logic;
	             x : in x_matrix_pwas;
	        scaleA : in coeff_matrix_pwas;
	            mu : in mu_matrix_pwas;
	             w : in mem_matrix_pwas;
	      mul_done : in std_logic;
	mul_done_scale : out std_logic;
	mul_done_final : out std_logic;
	        mul_x1 : out mul_in_matrix_pwas;
	        mul_x2 : out mul_in_matrix_pwas;
	     mul_start : out std_logic);
end mulBankManager;

architecture Behavioral of mulBankManager is

	signal state : integer range -1 to N_FUN_PWAS+1;

begin

	proc_state: process(clk,reset)
	begin
		if reset = '0' then
			state <= -1;
		elsif rising_edge(clk) then
			if state = -1 then
				if sorter_done = '1' then
					state <= state+1;
				else
					state <= state;
				end if;
			elsif state = 0 then
				if mem_done = '1' then
					state <= state+1;	
				else
					state <= state;
				end if;
			elsif state = N_FUN_PWAS+1 then
				state <= -1;
			else 
				state <= state+1;			
			end if;
		end if;
	end process;
	
	proc_mulBank: process(state,x,mu,w,scaleA)						
	begin
		if state = -1 then
			for i in 0 to N_DIM_PWAS-1 loop
				mul_x1(i) <= signed(resize(unsigned(x(i)),N_BIT_MUL_PWAS));
				mul_x2(i) <= resize(scaleA(i),N_BIT_MUL_PWAS);
			end loop;		
			mul_x1(N_DIM_PWAS) <= (others => '0');
			mul_x2(N_DIM_PWAS) <= (others => '0');
		else
			for i in 0 to N_DIM_PWAS loop
				mul_x1(i) <= resize("0"&signed(mu(i)(N_BIT_COEFF_PWAS-1 downto 1)),N_BIT_MUL_PWAS);
			end loop;	
			if state < N_FUN_PWAS then
				for i in 0 to N_DIM_PWAS loop
					mul_x2(i) <= resize(w(state)(i),N_BIT_MUL_PWAS);
				end loop;	
			else
				for i in 0 to N_DIM_PWAS loop
					mul_x2(i) <= resize(w(N_FUN_PWAS-1)(i),N_BIT_MUL_PWAS);
				end loop;	
			end if;
		end if;
	end process;
	
	mul_done_scale <= mul_done when state = -1 else
	                  '0';
	
	mul_done_final <= mul_done when state > -1 else
	                  '0';				
   
	mul_start <= start when state = -1 else
	             '1' when (state = 0 and mem_done = '1') else
	             '1' when (state > 0 and state < N_FUN_PWAS) else
	             '0';		

end Behavioral;

