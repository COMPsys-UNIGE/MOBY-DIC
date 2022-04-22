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

entity Scale is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	     start : in std_logic;
	         x : in x_matrix_pwas;
	      Axmb : in mul_out_matrix_pwas;
	       xmb : out x_matrix_pwas;
	    scaleA : out coeff_matrix_pwas;
	     x_int : out int_matrix_pwas;
	     x_dec : out dec_matrix_pwas;
	      done : out  std_logic);
end Scale;

architecture Behavioral of Scale is

begin

	proc_scaleB: process(reset,x)
	begin
		if reset = '0' then
			for i in 0 to N_DIM_PWAS-1 loop
				xmb(i) <= (others => '0');
			end loop;
		else
			for i in 0 to N_DIM_PWAS-1 loop
				xmb(i) <= x(i)-scaleInput_b_PWAS(i);
			end loop;
		end if;
	end process;


	proc_intdec: process(clk,start,reset)
		variable signal_Axmb : mul_out_matrix_pwas;
	begin
		if reset = '0' then
			for i in 0 to N_DIM_PWAS-1 loop
				x_int(i) <= (others => '0');
				x_dec(i) <= (others => '0');
			end loop;
		elsif rising_edge(clk) and start = '1' then
			for i in 0 to N_DIM_PWAS-1 loop
				if Axmb(i) > NP_PWAS(i) then
					signal_Axmb(i) := NP_PWAS(i);
				elsif Axmb(i) < ZERO_PWAS then
					signal_Axmb(i) := ZERO_PWAS;
				else
					signal_Axmb(i) := Axmb(i);
				end if;
				x_int(i) <= resize(unsigned(signal_Axmb(i)(POINT_POSITION_PWAS(i)+N_BIT_INT_PWAS(i)-1 downto POINT_POSITION_PWAS(i))),N_INT_MAX_PWAS);
				x_dec(i) <= unsigned(signal_Axmb(i)(POINT_POSITION_PWAS(i)-1 downto POINT_POSITION_PWAS(i)-N_BIT_COEFF_PWAS));
			end loop;
		end if;
	end process;

	scaleA <= scaleInput_A_PWAS;

	done <= '0' when reset = '0' else
	        start when rising_edge(clk);

end Behavioral;

