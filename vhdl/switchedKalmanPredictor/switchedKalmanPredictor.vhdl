----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    09:04:34 01/07/2016 
-- Design Name: 
-- Module Name:    kalmanPredictor - Behavioral 
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
use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity switchedKalmanPredictor is
	Port( clk : in std_logic;
			reset : in std_logic;
			start : in std_logic;
			stop : out std_logic;
			mulIn1 : out mul_in_matrix_OBS;
			mulIn2 : out mul_in_matrix_OBS;
			mulOut : in mul_out_matrix_OBS;
            dynamic : in integer range 0 to nDyn_OBS-1;
            in_matrix : in compositeInputType_OBS;
            out_matrix : out compositeStateType_OBS
);
end switchedKalmanPredictor;

architecture Behavioral of switchedKalmanPredictor is

component statePredict is
    Port ( clk : in  STD_LOGIC;
           reset : in  STD_LOGIC;
			  start : in std_logic;
			  valid_output : out std_logic;
			  mulIn1 : out mul_in_matrix_OBS;
			  mulIn2 : out mul_in_matrix_OBS;
			  mulOut : in mul_out_matrix_OBS;
              dynamic : in integer range 0 to nDyn_OBS-1;
           x : in  compositeStateType_OBS; -- x = [x;d]
           u : in  compositeInputType_OBS; -- u = [u;p;y]
			  x_next : out compositeStateType_OBS);
end component;

signal u_in : compositeInputType_OBS;

signal x_next_reg : compositeStateType_OBS;
signal x_next_tmp : compositeStateType_OBS;
signal valid_output_signal : std_logic;
signal stop_reg : std_logic;

begin

compute_port_map : statePredict port map(clk => clk,
						reset => reset,
						start => start,
						valid_output => valid_output_signal,
						mulIn1 => mulIn1,
						mulIn2 => mulIn2,
						mulOut => mulOut,
                        dynamic => dynamic,
						x => x_next_reg,
						u => u_in,
						x_next => x_next_tmp);

-- x_next reg process
x_next_reg_process : process(reset,clk)
begin 
if reset = '0' then 
	for i in 0 to nx_OBS+nd_OBS-1 loop
		x_next_reg(i) <= initialState(i);
	end loop;
elsif rising_edge(clk) then 
	if valid_output_signal = '1' then
		x_next_reg <= x_next_tmp;
	end if;
end if;
end process;

-- stop reg process; stop is active one clock cycle after enable
stop_reg_process : process(reset,clk)
begin 
if reset = '0' then 
	stop_reg <= '0';
elsif rising_edge(clk) then 
	stop_reg <= valid_output_signal;
end if;
end process;

out_matrix <= x_next_reg;
u_in <= in_matrix; 


stop <= stop_reg;

end Behavioral;

