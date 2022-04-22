----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    15:26:43 01/07/2016 
-- Design Name: 
-- Module Name:    kalmanFilter - Behavioral 
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

entity switchedKalmanFilter is
	Port( clk : in std_logic;
			reset : in std_logic;
			start_update : in std_logic;
			start_predict : in std_logic;
			stop_update : out std_logic;
			stop_predict : out std_logic;
            dynamic : in integer range 0 to nDyn_OBS-1;
			mulIn1 : out mul_in_matrix_OBS;
			mulIn2 : out mul_in_matrix_OBS;
			mulOut : in mul_out_matrix_OBS;
            in_matrix : in compositeInputType_OBS;
            out_matrix_pred : out compositeStateType_OBS;
            out_matrix_update : out compositeStateType_OBS
);
end switchedKalmanFilter;


architecture Behavioral of switchedKalmanFilter is

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

component stateUpdate is
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
			  x_update : out compositeStateType_OBS);
end component;

signal u_in : compositeInputType_OBS;

type FilterFSMstate is (Idle, Predict, Update);
signal current_state, next_state : FilterFSMstate;


signal x_next_reg : compositeStateType_OBS;
signal x_next_tmp : compositeStateType_OBS;
signal x_updated_reg : compositeStateType_OBS;
signal x_updated_tmp : compositeStateType_OBS;
signal valid_output_update_signal : std_logic;
signal valid_output_predict_signal : std_logic;
signal stop_update_reg : std_logic;
signal stop_predict_reg : std_logic;

--signal start_update_signal : std_logic;
--signal start_predict_signal : std_logic;

signal mux_sel : std_logic;

signal mulIn1_predict : mul_in_matrix_OBS;
signal mulIn2_predict : mul_in_matrix_OBS;


signal mulIn1_update : mul_in_matrix_OBS;
signal mulIn2_update : mul_in_matrix_OBS;


begin

---------- internal FSM process ------------------

state_proc : process(clk,reset)
begin
	if reset = '0' then
		current_state <= Idle;
	elsif rising_edge(clk) then
		current_state <= next_state;
	end if;
end process;

next_state_proc : process(current_state,valid_output_predict_signal,valid_output_update_signal,start_update,start_predict)
begin
case current_state is
  when Idle =>   
		if start_update = '1' then
			next_state <= Update;
		elsif start_predict = '1' then
			next_state <= Predict;
		else
			next_state <= Idle;
		end if;
  when Update =>   
		if valid_output_update_signal = '1' then
			next_state <= Idle;
		else 
			next_state <= update;
		end if;
  when Predict =>   
		if valid_output_predict_signal = '1' then
			next_state <= Idle;
		else 
			next_state <= predict;
		end if;
  when others => 
	next_state <= Idle;
end case;
end process; 

mux_sel_proc : process(current_state,start_update,start_predict)
begin
case current_state is
  when Idle =>   
		if start_update = '1' then
			mux_sel <= '0';
			elsif start_predict = '1' then
			mux_sel <= '1';
			else
			mux_sel <= '0';
		end if;
  when Predict =>   
		mux_sel <= '1';
  when others => 
	mux_sel <= '0';
end case;
end process;

--start_update_signal <= start_update when current_state = Idle 
--							else '0';
--start_predict_signal <= start_predict when current_state = Idle 
--							else '0';


---------- end of internal FSM process ------------------



-------- MUX multiplier access -------------------------
mulIn1 <= mulIn1_Predict when mux_sel = '1' 
			else mulIn1_Update;
			
mulIn2 <= mulIn2_Predict when mux_sel = '1' 
			else mulIn2_Update;
-------- end of MUX multiplier access -------------------

predict_port_map : statePredict port map(clk => clk,
						reset => reset,
						start => start_predict,
						valid_output => valid_output_predict_signal,
						mulIn1 => mulIn1_predict,
						mulIn2 => mulIn2_predict,
						mulOut => mulOut,
                        dynamic => dynamic,
						x => x_next_reg,
						u => u_in,
						x_next => x_next_tmp);
						
update_port_map : stateUpdate port map(clk => clk,
						reset => reset,
						start => start_update,
						valid_output => valid_output_update_signal,
						mulIn1 => mulIn1_update,
						mulIn2 => mulIn2_update,
						mulOut => mulOut,
                        dynamic => dynamic,
						x => x_next_reg,
						u => u_in,
						x_update => x_updated_tmp);

-- x_next reg process
x_next_reg_process : process(reset,clk)
begin 
if reset = '0' then 
	for i in 0 to nx_OBS+nd_OBS-1 loop
		x_next_reg(i) <= initialState(i);
	end loop;
elsif rising_edge(clk) then 
	if valid_output_predict_signal = '1' then
		x_next_reg <= x_next_tmp;
	end if;
end if;
end process;


-- x_update reg process
x_update_reg_process : process(reset,clk)
begin 
if reset = '0' then 
	for i in 0 to nx_OBS+nd_OBS-1 loop
		x_updated_reg(i) <= initialState(i);
	end loop;
elsif rising_edge(clk) then 
	if valid_output_update_signal = '1' then
		x_updated_reg <= x_updated_tmp;
	end if;
end if;
end process;

-- stop reg process; stop is active one clock cycle after enable
stop_reg_process : process(reset,clk)
begin 
if reset = '0' then 
	stop_update_reg <= '0';
	stop_predict_reg <= '0';
elsif rising_edge(clk) then 
	stop_update_reg <= valid_output_update_signal;
	stop_predict_reg <= valid_output_predict_signal;
end if;
end process;

out_matrix_update <= x_updated_reg;
out_matrix_pred <= x_next_reg;
u_in <= in_matrix;

stop_update <= stop_update_reg;
stop_predict <= stop_predict_reg;

end Behavioral;