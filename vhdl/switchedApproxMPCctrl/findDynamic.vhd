----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    14:47:01 02/01/2016 
-- Design Name: 
-- Module Name:    findDynamic - Behavioral 
-- Project Name: 
-- Target Devices: 
-- Tool versions: 
-- Description: 
--
-- Dependencies: 
--
-- Revision: 
-- Revision 0.01 - File Created
-- Additional Comments: 
--
----------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
use IEEE.NUMERIC_STD.ALL;
use work.controller_package.all;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity findDynamic is
    Port ( clk : in  STD_LOGIC;
           reset : in  STD_LOGIC;
			  start : in std_logic;
			  valid_output : out std_logic;
			  mulIn1 : out mul_in_matrix_CTRL;
			  mulIn2 : out mul_in_matrix_CTRL;
			  mulOut : in mul_out_matrix_CTRL;
			  x : in  x_matrix_CTRL;
           dynamic : out integer range 0 to N_DYN_CTRL-1);
end findDynamic;

architecture Behavioral of findDynamic is

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

signal counter : integer range 0 to N_EDGE+2;

signal dynamicOut : integer range 0 to N_DYN_CTRL-1;
signal dynamicOut_next : integer range 0 to N_DYN_CTRL-1;

signal mul_out_reg : mul_out_matrix_CTRL;
signal result_reg : std_logic_vector(N_EDGE-1 downto 0); -- vector containing result of inequalities
signal Kreg : signed(N_BIT_COEFF_CTRL-1 downto 0);

begin

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

dynamicOut_reg_process : process(clk,reset)
begin
if reset = '0' then
	dynamicOut <= 0;
elsif rising_edge(clk) then
if counter = N_EDGE+1 then
	dynamicOut <= dynamicOut_next;
end if;
end if;
end process;

sum_reg_process : process(clk,reset)
variable tmp : signed(maxSchift_findDyn+N_BIT_COEFF_CTRL+N_BIT_CTRL+N_DIM_CTRL-1 downto 0);
begin
if reset = '0' then
	result_reg <= (others => '0');
elsif rising_edge(clk) then
tmp := (others => '0');
	for i in 0 to N_DIM_CTRL-1 loop
		tmp := tmp + SHIFT_LEFT(resize(mul_out_reg(i),maxSchift_findDyn+N_BIT_COEFF_CTRL+N_BIT_CTRL+N_DIM_CTRL),shift_findDyn(i));
	end loop;
tmp := tmp + SHIFT_LEFT(resize(Kreg,maxSchift_findDyn+N_BIT_COEFF_CTRL+N_BIT_CTRL+N_DIM_CTRL),shift_findDyn(N_DIM_CTRL));

if counter > 0 and counter < N_EDGE+1 then
	result_reg(counter-1) <= std_logic(tmp(maxSchift_findDyn+N_BIT_COEFF_CTRL+N_BIT_CTRL+N_DIM_CTRL-1));
end if;

end if;
end process;

Kreg_process : process(clk,reset)
begin
if reset = '0' then
	Kreg <= (others => '0');
elsif rising_edge(clk) then
	if counter < N_EDGE then
		Kreg <= HK(counter)(N_DIM_CTRL);
	else
		Kreg <= HK(N_EDGE-1)(N_DIM_CTRL);
	end if;
end if;
end process;


mul_Out_reg_process : process(clk,reset)
begin
if reset = '0' then
	for i in 0 to N_DIM_CTRL-1 loop
		mul_out_reg(i) <= (others => '0');
	end loop;
elsif rising_edge(clk) then
	mul_out_reg <= mulOut;
end if;
end process;

mul1_in_process : process(x)
begin
for i in 0 to N_DIM_CTRL-1 loop
	mulIn1(i) <= resize(signed(x(i)),N_BIT_MUL_CTRL);
end loop;
end process;

mul2_in_process : process(counter)
begin
if counter < N_EDGE then
for i in 0 to N_DIM_CTRL-1 loop
	mulIn2(i) <= resize(HK(counter)(i),N_BIT_MUL_CTRL);
end loop;
else
for i in 0 to N_DIM_CTRL-1 loop
	mulIn2(i) <= resize(HK(N_EDGE-1)(i),N_BIT_MUL_CTRL);
end loop;
end if;
end process;

counter_process : process(clk,reset)
begin
if reset = '0' then
	counter <= 0;
elsif rising_edge(clk) then
	if counter = 0 then
		if start = '1' then
			counter <= counter + 1;
		end if;
	elsif counter = N_EDGE+2 then
		counter <= 0;
	else
			counter <= counter + 1;
	end if;
end if;
end process;



dynamic <= dynamicOut;

valid_output <= '1' when counter = N_EDGE+2
					else '0';

end Behavioral;

