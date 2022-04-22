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
use work.pwagFunctionPackage.all;

entity compute is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	      leaf : in std_logic;
	   restart : in std_logic;
		   addr : in addr_array_PWAG;
	     mul_y : in mul_out_matrix_PWAG;
	    mul_x2 : out mul_in_matrix_PWAG;
	         y : out y_matrix_PWAG;
      decision: out std_logic;
 compute_done : out  std_logic;
 compute_func : out std_logic);
end compute;

-- Architecture declaration

architecture Behavioral of compute is

component Memory is
	Port ( addr : in unsigned(N_BIT_ADDR_PWAG-1 downto 0);
		clk : in std_logic;
		res : out memory_out_matrix_PWAG);
end component;

signal compute_func_reg : std_logic;
signal compute_func_next : std_logic;

signal function_counter : integer range 0 to N_FUN_PWAG+3;
signal calc_counter : integer range 0 to 2;

signal addr_from_FSM : addr_array_PWAG;
signal addr_to_memory : unsigned(N_BIT_ADDR_PWAG-1 downto 0);

signal mul_y_reg : mul_out_matrix_PWAG;

signal Ax_plus_B_tmp : signed(N_BIT_ADDER_PWAG-1 downto 0);
signal Ax_plus_B_reg : signed(N_BIT_ADDER_PWAG-1 downto 0);

signal rescale_sum : signed(N_BIT_COEFF_PWAG+N_BIT_ADDER_PWAG-1 downto 0);

signal out_reg : y_matrix_PWAG;


signal memory_out : memory_out_matrix_PWAG;

signal memReg : signed(N_BIT_COEFF_PWAG-1 downto 0);

begin	

memory_portmap : Memory 
	port map( addr => addr_to_memory,
		clk => clk,
		res => memory_out);


calc_counter_process : process(clk,reset)
begin
if reset = '0' then
		calc_counter <= 0;
elsif rising_edge(clk) then
	if restart = '1' or leaf = '1' then
		calc_counter <= 0;
	else
		if calc_counter = 2 then
			calc_counter <= 0;
		else
			calc_counter <= calc_counter + 1;
		end if;
	end if;
end if;
end process;

function_counter_process : process(clk,reset)
begin
if reset = '0' then
		function_counter <= 0;
elsif rising_edge(clk) then
if restart = '1' then
function_counter <= 0;
else
	if function_counter = N_FUN_PWAG+3 then
		function_counter <= N_FUN_PWAG+3;
	else
		if leaf = '1' then
				function_counter <= function_counter + 1;	
		else
			function_counter <= 0;
		end if;
	end if;
end if;
end if;
end process;
	
y_mul_reg_process : process(reset,clk)
begin
if reset = '0' then
for i in 0 to N_DIM_PWAG-1 loop
	mul_y_reg(i) <= (others => '0');
end loop;
elsif rising_edge(clk) then
	mul_y_reg <= mul_y;
end if;
end process;	

mul_in2_signal_process : process(memory_out)
begin
for i in 0 to N_DIM_PWAG-1 loop
	mul_x2(i) <= resize(memory_out(i),N_BIT_MUL_PWAG);
end loop;
end process;

sum_process : process(mul_y_reg,memReg)
variable tmp : signed(N_BIT_ADDER_PWAG-1 downto 0);
begin
tmp := (others => '0');
for i in 0 to N_DIM_PWAG-1 loop
	tmp := tmp + SHIFT_LEFT(resize(mul_y_reg(i)(N_BIT_PWAG+N_BIT_COEFF_PWAG-1 downto 0),N_BIT_ADDER_PWAG),preshift(i));
end loop;
tmp := tmp + SHIFT_LEFT(resize(memReg,N_BIT_ADDER_PWAG),preshift(N_DIM_PWAG));
Ax_plus_B_tmp <= tmp;
end process;

constant_mem_reg : process(clk,reset)
begin
if reset = '0' then
	memReg <= (others => '0');
elsif rising_edge(clk) then
	memReg <= 	memory_out(N_DIM_PWAG);
end if;
end process;


sum_reg_process : process(clk,reset)
begin
if reset = '0' then
		Ax_plus_B_reg <= (others => '0');
elsif rising_edge(clk) then
	Ax_plus_B_reg <= Ax_plus_B_tmp;
end if;
end process;


reg_out_process : process(clk,reset)
begin
if reset = '0' then
	for i in 0 to N_FUN_PWAG-1 loop
		out_reg(i) <= (others => '0');
	end loop;
elsif rising_edge(clk) then
	if function_counter >= 3 and function_counter < N_FUN_PWAG+3 then
	--saturation
	if rescale_sum > maxVal_PWAG then
		out_reg(function_counter-3) <= maxVal_PWAG(rescaleOutDecPosistion_PWAG+N_BIT_OUT_PWAG-1 downto rescaleOutDecPosistion_PWAG);
	elsif rescale_sum < minVal_PWAG then
		out_reg(function_counter-3) <= minVal_PWAG(rescaleOutDecPosistion_PWAG+N_BIT_OUT_PWAG-1 downto rescaleOutDecPosistion_PWAG);
	else
	out_reg(function_counter-3) <= rescale_sum(rescaleOutDecPosistion_PWAG+N_BIT_OUT_PWAG-1 downto rescaleOutDecPosistion_PWAG);--out_reg;
	end if;
	end if;
end if;
end process;

compute_func_process : process(clk,reset)
begin
if reset = '0' then
	compute_func_reg <= '0';
elsif rising_edge(clk) then
	compute_func_reg <= compute_func_next;
end if;
end process;

compute_func_next <= '1' when function_counter = N_FUN_PWAG+2
else '0';

compute_func <= compute_func_reg;

decision <= Ax_plus_B_tmp(N_BIT_ADDER_PWAG-1);

rescale_sum <= Ax_plus_B_reg*alpha_PWAG+beta_PWAG;

compute_done <= '1' when 	calc_counter = 2
		else '0';
	
addr_from_FSM <= addr;

addr_to_memory <= addr_from_FSM(function_counter) when function_counter < N_FUN_PWAG
	else addr_from_FSM(N_FUN_PWAG-1);

y <= out_reg;

end Behavioral;