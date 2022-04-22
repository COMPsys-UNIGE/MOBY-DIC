----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    08:37:15 01/12/2016 
-- Design Name: 
-- Module Name:    kalmanObserver - Behavioral 
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

use work.observer_package.all;
-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity kalmanPredictorInterface is
Port(clk : in std_logic;
      ce : in std_logic;
	reset : in std_logic;
    sample : out std_logic;
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---
		);
end kalmanPredictorInterface;

architecture Behavioral of kalmanPredictorInterface is

component kalmanPredictor is
	Port( clk : in std_logic;
			reset : in std_logic;
			start : in std_logic;
			stop : out std_logic;
			  mulIn1 : out mul_in_matrix_OBS;
			  mulIn2 : out mul_in_matrix_OBS;
			  mulOut : in mul_out_matrix_OBS;
            in_matrix : in compositeInputType_OBS;
            out_matrix : out compositeStateType_OBS
);
end component;

component mulBank is
Port(			  x1 : in mul_in_matrix_OBS;
			  x2 : in mul_in_matrix_OBS;
			  y : out mul_out_matrix_OBS);
end component;

component inputShift is
Port(
u_in : in compositeInputType_OBS;
	u_out : out compositeInputType_OBS);
end component;

component outputShift is
Port(
x_in : in compositeStateType_OBS;
	x_out : out compositeStateType_OBS);
end component;

signal counter : integer range 0 to sampling_latency_OBS-1;
signal start : std_logic;
signal stop : std_logic;

signal mulIn1 : mul_in_matrix_OBS;
signal mulIn2 : mul_in_matrix_OBS;
signal mulOut : mul_out_matrix_OBS; 

signal	x_out_reg : compositeStateType_OBS;
signal	x_out : compositeStateType_OBS;
signal	x_out_scaled : compositeStateType_OBS;
signal u_in_reg : compositeInputType_OBS;
signal u_in : compositeInputType_OBS;
signal u_in_scaled : compositeInputType_OBS;

begin

kalmanPredictor_portmap : kalmanPredictor port map
	(clk => clk,
	reset => reset,
	start => start,
	mulOut => mulOut,
	mulIn1 => mulIn1,
	mulIn2 => mulIn2,
    in_matrix => u_in_reg,
    out_matrix => x_out);

inputShift_portmap : inputShift port map
    (u_in => u_in,
    u_out => u_in_scaled);

outputShift_portmap : outputShift port map
    (x_in => x_out,
    x_out => x_out_scaled);


mulBank_portmap : mulBank port map
	(	x1 => mulIn1,
	x2 => mulIn2,
	y => mulOut);

counter_process : process (clk,reset)
begin
if reset = '0' then
	counter <= 0;
elsif rising_edge(clk) then
	if counter = sampling_latency_OBS-1 then 
		counter <= 0;
	else
		counter <= counter + 1;
	end if;
end if;
end process;

start <= '1' when counter = 1
		else '0';

sample <= '1' when counter = 0
		else '0';

in_reg : process(reset,clk)
begin
if reset = '0' then
for i in 0 to nu_OBS+np_OBS+ny_OBS-1 loop
    u_in_reg(i) <= (others => '0');
end loop;
elsif rising_edge(clk) then
if counter = 0 then
    u_in_reg <= u_in_scaled;
end if;
end if;
end process;

out_reg : process(reset,clk)
begin
if reset = '0' then
for i in 0 to nx_OBS+nd_OBS-1 loop
    x_out_reg(i) <= (others => '0');
end loop;
elsif rising_edge(clk) then
if counter = 0 then
    x_out_reg <= x_out_scaled;
end if;
end if;
end process;

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

end Behavioral;

