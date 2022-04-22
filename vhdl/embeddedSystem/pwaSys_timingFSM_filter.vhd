----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    11:12:17 01/26/2016 
-- Design Name: 
-- Module Name:    timingFSM - Behavioral 
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

use work.embeddedSystemPackage.all;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity timingFSM is
Port(clk : in std_logic;
	reset : in std_logic;
	sample : out std_logic;
    update : out std_logic;
	end_update : in std_logic;
	ctrl : out std_logic;
	end_control : in std_logic;
	predict : out std_logic;
	end_predict : in std_logic;
    findDyn : out std_logic;
    end_findDyn : in std_logic;
	mux_select : out std_logic_vector(1 downto 0);
	sel_x_in_dyn : out std_logic
	);
end timingFSM;

architecture Behavioral of timingFSM is

type FSMstate is (Idle,findUpdateDynStart, findUpdateDynState, UpdateStart,UpdateState, findCtrlDynStart, findCtrlDynState,  ControlStart, ControlState, PredictStart,PredictState);

signal state,nextState : FSMstate;

signal obs_counter : integer range 0 to observer_latency_ES-1;
signal ctrl_counter : integer range 0 to controller_mul_observer_latency_ES-1;

begin

counter_proc : process(clk,reset)
begin
if reset = '0' then
obs_counter <= 0;
ctrl_counter <= 0;
elsif rising_edge(clk) then
	if obs_counter = observer_latency_ES-1 then
		obs_counter <= 0;
		if ctrl_counter = controller_mul_observer_latency_ES-1 then
			ctrl_counter <= 0;
		else
			ctrl_counter <= ctrl_counter + 1;
		end if;
	else
		obs_counter <= obs_counter + 1;
	end if;
end if;
end process;

state_reg : process(clk,reset)
begin
if reset = '0' then
	state <= Idle;
elsif rising_edge(clk) then
	state <= nextState;
end if;
end process;

next_state_proc : process(state,obs_counter,ctrl_counter,end_predict,end_control,end_update,end_findDyn)
begin
case state is
	when Idle => 
		if obs_counter = 0 then
            if ctrl_counter = 0 then
                nextState <= findUpdateDynStart;
            else
                nextState <= PredictStart;
            end if;
		else
			nextState <= Idle;
		end if;
    when findUpdateDynStart => 
        nextState <= findUpdateDynState;
    when findUpdateDynState =>
        if end_findDyn = '1' then
            nextState <= UpdateStart;
        else
            nextState <= findUpdateDynState;
        end if;
    when UpdateStart => 
        nextState <= UpdateState;
    when UpdateState =>
        if end_update = '1' then
            nextState <= findCtrlDynStart;
        else
            nextState <= UpdateState;
        end if;
    when findCtrlDynStart => 
        nextState <= findCtrlDynState;
    when findCtrlDynState =>
        if end_findDyn = '1' then
			if ctrl_counter = 0 then
				nextState <= ControlStart;
			else
				nextState <= PredictStart;
			end if;
		else
			nextState <= findCtrlDynState;
		end if;
	when ControlStart =>
		nextState <= ControlState;
	when ControlState =>
		if end_control = '1' then
			nextState <= PredictStart;
		else
			nextState <= ControlState;
		end if;
	when PredictStart =>
		nextState <= PredictState;
	when PredictState =>
		if end_predict = '1' then 
			nextState <= Idle;
		else
			nextState <= PredictState;
		end if;
	when others =>
		nextState <= Idle;
end case;
end process;

output_process : process(state)
begin
ctrl <= '0';
predict <= '0';
update <= '0';
findDyn <= '0';
mux_select <= "00";
sel_x_in_dyn <= '0';
case state is
	when ControlStart =>
			ctrl <= '1';
			mux_select <= "01";
	when ControlState =>
			mux_select <= "01";
	when PredictStart =>
			predict <= '1';
    when UpdateStart =>
			update <= '1';
    when findCtrlDynStart => 
            findDyn <= '1';
            mux_select <= "10";
    when findCtrlDynState => 
            mux_select <= "10";
    when findUpdateDynStart => 
            findDyn <= '1';
            mux_select <= "10";
				sel_x_in_dyn <= '1';
    when findUpdateDynState => 
            mux_select <= "10";
				sel_x_in_dyn <= '1';
	when others =>
end case;
end process;

sample <= '1' when obs_counter = observer_latency_ES-1 
	else '0';

end Behavioral;

