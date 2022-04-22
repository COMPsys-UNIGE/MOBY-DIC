----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    09:04:34 01/07/2016 
-- Design Name: 
-- Module Name:    controllerCircuit - Behavioral 
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


use work.controller_package.all;
-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

architecture Behavioral of switchedControllerInterface is

component findDynamic is
    Port ( clk : in  STD_LOGIC;
           reset : in  STD_LOGIC;
			  start : in std_logic;
			  valid_output : out std_logic;
			  mulIn1 : out mul_in_matrix_CTRL;
			  mulIn2 : out mul_in_matrix_CTRL;
			  mulOut : in mul_out_matrix_CTRL;
			  x : in  x_matrix_CTRL;
           dynamic : out integer range 0 to N_DYN_CTRL-1);
end component;


component switchedControllerCircuit is
	Port( clk : in std_logic;
			reset : in std_logic;
			start : in std_logic;
			 done : out std_logic;
          dynamic : in integer range 0 to N_DYN_CTRL-1;
                x : in x_matrix_CTRL;
            mul_y : in mul_out_matrix_CTRL;
           mul_x1 : out mul_in_matrix_CTRL;
           mul_x2 : out mul_in_matrix_CTRL;
                u : out u_matrix_CTRL
);

end component;

component mulBank is
Port(			  x1 : in mul_in_matrix_CTRL;
			  x2 : in mul_in_matrix_CTRL;
			  y : out mul_out_matrix_CTRL);
end component;

	
	-- Signals declaration
	
	signal signal_u: u_matrix_CTRL;
	signal signal_x: x_matrix_CTRL;


signal counter : integer range 0 to sampling_latency_CTRL-1;
signal u_reg : u_matrix_CTRL;
signal x_reg : x_matrix_CTRL;


signal start_pwas : std_logic;
signal start_findDynamic : std_logic;
signal sign_sample : std_logic;
signal signal_done : std_logic;
 
signal mulIn1 : mul_in_matrix_CTRL;
signal mulIn2 : mul_in_matrix_CTRL;
signal mulOut : mul_out_matrix_CTRL;

signal mulIn1_pwas : mul_in_matrix_CTRL;
signal mulIn2_pwas : mul_in_matrix_CTRL;

signal mulIn1_dyn : mul_in_matrix_CTRL;
signal mulIn2_dyn : mul_in_matrix_CTRL;

signal  dynamic_signal : integer range 0 to N_DYN_CTRL-1;
signal mux_sel : std_logic;
signal valid_dyn : std_logic;
begin

	-- Port maps
	
	inst_controllerCircuit: switchedControllerCircuit
		port map (  clk => clk,
		         reset => reset,
		         start => start_pwas,
                dynamic => dynamic_signal,
		             x => x_reg,
		         mul_y => mulOut,
		        mul_x1 => mulIn1_pwas,
		        mul_x2 => mulIn2_pwas,
		             u => signal_u,
		          done => signal_done);

inst_findDynamic : findDynamic 
    Port map ( clk => clk,
           reset => reset,
			  start => start_findDynamic,
			  valid_output => valid_dyn,
			  mulIn1 => mulIn1_dyn,
			  mulIn2 => mulIn2_dyn,
			  mulOut => mulOut,
			  x => x_reg,
           dynamic => dynamic_signal);

inst_mulBank : mulBank 
    port map(x1 => mulIn1,
			  x2 => mulIn2,
			  y => mulOut);


counter_process : process(clk,reset)
begin
if reset = '0' then
    counter <= 0;
elsif rising_edge(clk) then
    if counter = sampling_latency_CTRL-1 then
        counter <= 0;
    else
        counter <= counter+1;
    end if;
end if;
end process;



sign_sample <= '1' when counter = 0
    else '0';

start_findDynamic <= '1' when counter = 1
    else '0';

input_reg_process : process(clk,reset)
begin
if reset = '0' then
    for i in 0 to nx_CTRL+npar_CTRL+nd_CTRL+nxref_CTRL-1 loop
        x_reg(i) <= (others => '0');
    end loop;
elsif rising_edge(clk) then
    if sign_sample = '1' then
        x_reg <= signal_x;
    end if;
end if;
end process;


output_reg_process : process(clk,reset)
begin
if reset = '0' then
    for i in 0 to nu_CTRL-1 loop
        u_reg(i) <= defaultOutput(i);
    end loop;
elsif rising_edge(clk) then
    if signal_done = '1' then
        u_reg <= signal_u;
    end if;
end if;
end process;


mul_process : process(mulIn1_pwas,mulIn2_pwas,mulIn1_dyn,mulIn2_dyn,mux_sel)
begin
if mux_sel = '1' then
        mulIn1 <= mulIn1_pwas;
        mulIn2 <= mulIn2_pwas;
else
        mulIn1 <= mulIn1_dyn;
        mulIn2 <= mulIn2_dyn;
end if;
end process;


mux_sel_process : process(clk,reset)
begin
if reset = '0' then
    mux_sel <= '0';
    start_pwas <= '0';
elsif rising_edge(clk) then
    start_pwas <= '0';
    if valid_dyn = '1' then
        mux_sel <= '1';
        start_pwas <= '1';
    elsif signal_done = '1' then
        mux_sel <= '0';
    end if;
end if;
end process;

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

sample <= sign_sample;
done <= signal_done;
end Behavioral;

