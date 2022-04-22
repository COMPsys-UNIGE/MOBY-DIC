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

entity switchedKalmanFilterInterface is
Port(clk : in std_logic;
      ce : in std_logic;
	reset : in std_logic;
    sample : out std_logic;
    u_ready : in std_logic;
    data_ready : out std_logic;
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---
		);
end switchedKalmanFilterInterface;

architecture Behavioral of switchedKalmanFilterInterface is

component switchedKalmanFilter is
	Port( clk : in std_logic;
			reset : in std_logic;
			start_update : in std_logic;
			start_predict : in std_logic;
			stop_update : out std_logic;
			stop_predict : out std_logic;
			mulIn1 : out mul_in_matrix_OBS;
			mulIn2 : out mul_in_matrix_OBS;
			mulOut : in mul_out_matrix_OBS;
			  dynamic : in integer range 0 to nDyn_OBS-1;
            in_matrix : in compositeInputType_OBS;
            out_matrix_pred : out compositeStateType_OBS;
            out_matrix_update : out compositeStateType_OBS
);
end component;

component mulBank is
Port(			  x1 : in mul_in_bigmatrix_OBS;
			  x2 : in mul_in_bigmatrix_OBS;
			  y : out mul_out_bigmatrix_OBS);
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

component findDynamic is
    Port ( clk : in  STD_LOGIC;
           reset : in  STD_LOGIC;
			  start : in std_logic;
			  valid_output : out std_logic;
			  mulIn1 : out mul_in_bigmatrix_OBS;
			  mulIn2 : out mul_in_bigmatrix_OBS;
			  mulOut : in mul_out_bigmatrix_OBS;
			  u : in  compositeInputType_OBS; -- u = [u;p;y]
			  x : in  compositeStateType_OBS;
           dynamic : out integer range 0 to nDyn_OBS-1);
end component;

signal signal_to_findDyn : compositeStateType_OBS; 

signal mulIn1_OBS : mul_in_matrix_OBS;
signal mulIn2_OBS : mul_in_matrix_OBS;
signal mulOut_OBS : mul_out_matrix_OBS;

signal mulIn1_dyn : mul_in_bigmatrix_OBS;
signal mulIn2_dyn : mul_in_bigmatrix_OBS;
signal mulOut_dyn : mul_out_bigmatrix_OBS;

signal mulIn1 : mul_in_bigmatrix_OBS;
signal mulIn2 : mul_in_bigmatrix_OBS;
signal mulOut : mul_out_bigmatrix_OBS;
 
signal dynamic_signal : integer range 0 to nDyn_OBS-1;

signal counter : integer range 0 to sampling_latency_OBS-1;
signal start_update : std_logic;
signal start_predict : std_logic;
signal stop_update : std_logic;
signal stop_predict : std_logic;
signal startFindDyn : std_logic;
signal stopFindDyn : std_logic;

signal	x_out_reg : compositeStateType_OBS;
signal	x_out : compositeStateType_OBS;
signal	x_out_scaled : compositeStateType_OBS;
signal	x_next_reg : compositeStateType_OBS;
signal	x_next : compositeStateType_OBS;
signal	x_next_scaled : compositeStateType_OBS;
signal u_in_reg : compositeInputType_OBS;
signal u_in : compositeInputType_OBS;
signal u_in_scaled : compositeInputType_OBS;

signal mux_Sel : std_logic; -- 0 findDynamic, 1 swObserver
signal toDynSel : std_logic;
type stateType is (Idle,startFindUpdate, findUpdate, startUpdate, Update, waitInput, startFindPredict, findPredict, startPredict, Predict);
signal state,next_state : stateType;


begin

findDynamic_portmap : findDynamic 
    Port map( clk => clk,
           reset => reset,
			  start => startFindDyn,
			  valid_output => stopFindDyn,
			  mulIn1 => mulIn1_dyn,
			  mulIn2 => mulIn2_dyn,
			  mulOut => mulOut_dyn,
			  u => u_in_reg,
			  x => signal_to_findDyn,
           dynamic => dynamic_signal);

signal_to_findDyn <= x_out_reg when toDynSel = '1'
	else x_next_reg;

kalmanFilter_portmap : switchedKalmanFilter port map
	(clk => clk,
	reset => reset,
	start_update => start_update,
    stop_update => stop_update,
    start_predict => start_predict,
	 stop_predict => stop_predict,
	mulIn1 => mulIn1_OBS,
	mulIn2 => 	 mulIn2_OBS,
	mulOut => mulOut_OBS,
	dynamic => dynamic_signal,
    in_matrix => u_in_reg,
    out_matrix_pred => x_next,
    out_matrix_update =>x_out);

inputShift_portmap : inputShift port map
    (u_in => u_in,
    u_out => u_in_scaled);

outputShiftUpdate_portmap : outputShift port map
    (x_in => x_out,
    x_out => x_out_scaled);
	 
outputShiftPred_portmap : outputShift port map
    (x_in => x_next,
    x_out => x_next_scaled);

mulBank_portmap : mulBank port map
	(	x1 => mulIn1,
	x2 => mulIn2,
	y => mulOut);


mul_process : process(mulOut,mulIn1_OBS,mulIn2_OBS,mulIn1_dyn,mulIn2_dyn,mux_sel)
begin
if mux_sel = '1' then
    for i in 0 to nd_OBS+nx_OBS-1 loop
        mulIn1(i) <= mulIn1_OBS(i);
        mulIn2(i) <= mulIn2_OBS(i);
    end loop;
else
    for i in 0 to nd_OBS+nx_OBS+np_OBS-1 loop
        mulIn1(i) <= mulIn1_dyn(i);
        mulIn2(i) <= mulIn2_dyn(i);
    end loop;
end if;

for i in 0 to nd_OBS+nx_OBS-1 loop
    mulOut_OBS(i) <= mulOut(i);
end loop;

for i in 0 to nd_OBS+nx_OBS+np_OBS-1 loop
    mulOut_Dyn(i) <= mulOut(i);
end loop;

end process;

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

sample <= '1' when counter = 0 
		else '0';

data_ready <= stop_update;


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
    if stop_update = '1' then
        x_out_reg <= x_out_scaled;
    end if;
end if;
end process;

pred_reg : process(reset,clk)
begin
if reset = '0' then
for i in 0 to nu_OBS+np_OBS+ny_OBS-1 loop
    x_next_reg(i) <= (others => '0');
end loop;
elsif rising_edge(clk) then
    if stop_predict = '1' then
        x_next_reg <= x_next_scaled;
    end if;
end if;
end process;


FSM_state_reg : process (clk,reset)
begin
if reset = '0' then
	state <= Idle;
elsif rising_edge(clk) then
	state <= next_state;
end if;
end process;

next_State_proc : process(state,counter,stopFindDyn,stop_update,stop_predict,u_ready)
begin
case state is
	when Idle =>
		if counter = 1 then
			next_state <= startFindUpdate;
		else
			next_state <= Idle;
		end if;
	when startFindUpdate =>
		next_state <= findUpdate;
	when findUpdate =>
		if stopFindDyn = '1' then
			next_state <= startUpdate;
		else
			next_state <= findUpdate;
		end if;
	when startUpdate =>
		next_state <= Update;
	when Update =>
		if stop_update = '1' then
			next_state <= waitInput;
		else
			next_state <= Update;
		end if;
	when waitInput  =>
		if u_ready = '1' then
			next_state <= startFindPredict;
		else
			next_state <= waitInput;
		end if;
	when startFindPredict =>
		next_state <= findPredict;
	when findPredict =>
		if stopFindDyn = '1' then
			next_state <= startPredict;
		else
			next_state <= findPredict;
		end if;
	when startPredict =>
		next_state <= Predict;
	when Predict =>
		if stop_predict = '1' then
			next_state <= Idle;
		else
			next_state <= Predict;
		end if;
	when others =>
end case;
end process;

output_FSM_proc : process(state)
begin
toDynSel <= '0';
mux_Sel <= '0';
startFindDyn <= '0';
start_update <= '0';
start_predict <= '0';
case state is
	when Idle =>
	when startFindUpdate =>
		startFindDyn <= '1';
		toDynSel <= '1';
	when findUpdate =>
		toDynSel <= '1';
	when startUpdate =>
		start_update <= '1';
		mux_Sel <= '1';
	when Update =>
		mux_Sel <= '1';
	when startFindPredict =>
		startFindDyn <= '1';
	when startPredict =>
		start_predict <= '1';
		mux_Sel <= '1';
	when Predict =>
		mux_Sel <= '1';
	when others =>
end case;
end process;

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

end Behavioral;

