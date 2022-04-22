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

entity mu_generator is
	Port (   clk : in std_logic;
	       reset : in std_logic;
	       start : in std_logic;
	x_dec_sorted : in  dec_matrix_pwas;
	          mu : out mu_matrix_pwas;
	        done : out std_logic);
end mu_generator;

-- Architecture declaration

architecture Behavioral of mu_generator is

	-- Constant = 1
	constant one : unsigned(N_BIT_COEFF_PWAS-1 downto 0) := (others => '1');
	
begin
	
	-- The weights are calculated by differences between the sorted fractional parts
	--
	-- mu(0) = 1-z_dec_sorted(0)
	-- ...
 	-- mu(i) = z_dec_sorted(i-1)-z_dec_sorted(i)
	-- ...
	-- mu(N_DIM_PWAS) = z_dec_sorted(N_DIM_PWAS-1)
	process(clk,reset,start,x_dec_sorted)
		variable mu_ext : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	begin
		if reset = '0' then
			for i in 0 to N_DIM_PWAS loop
				mu(i) <= (others => '0');
			end loop;
		elsif rising_edge(clk) and start = '1' then
			for i in 0 to N_DIM_PWAS loop
				if i = 0 then
					mu_ext := one-x_dec_sorted(0);
				elsif i = N_DIM_PWAS then
					mu_ext := x_dec_sorted(N_DIM_PWAS-1);
				elsif i < N_DIM_PWAS and i > 0 then
					mu_ext := x_dec_sorted(i-1)-x_dec_sorted(i);
				else
					mu_ext := (others=>'0');
				end if;
				mu(i) <= mu_ext(N_BIT_COEFF_PWAS-1 downto 0);
			end loop;
		end if;
	end process;
	
	done <= '0' when reset = '0' else
	        start when rising_edge(clk);

end Behavioral;
