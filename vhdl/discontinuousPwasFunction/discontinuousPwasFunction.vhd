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
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use work.discontinuousPwasFunctionPackage.all;

-- Entity declaration

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

-- Architecture declaration

architecture Behavioral of discontinuousPwasFunction is

	-- Components declaration
	
	component discontinuousPwasFunctionCircuit is
		Port ( clk : in std_logic;
		     reset : in std_logic;
		     start : in std_logic;
		         x : in x_matrix_pwas;
       dynamic : in integer range 0 to N_DYN_PWAS-1;
		     mul_y : in mul_out_matrix_pwas;
		    mul_x1 : out mul_in_matrix_pwas;
		    mul_x2 : out mul_in_matrix_pwas;
		         y : out y_matrix_pwas;
		      done : out std_logic);
	end component;

    component findDynamic is
Port ( clk : in  STD_LOGIC;
           reset : in  STD_LOGIC;
			  start : in std_logic;
			  valid_output : out std_logic;
			  mulIn1 : out mul_in_matrix_pwas;
			  mulIn2 : out mul_in_matrix_pwas;
			  mulOut : in mul_out_matrix_pwas;
			  x : in  x_matrix_pwas;
           dynamic : out integer range 0 to N_DYN_PWAS-1);
    end component;
	
	component mulBank is
		Port (  x1 : in mul_in_matrix_pwas;
		        x2 : in mul_in_matrix_pwas;
		         y : out mul_out_matrix_pwas);
	end component;
	
	-- Signals declaration
	
	signal signal_y: y_matrix_pwas;
	signal signal_x: x_matrix_pwas;
	signal signal_mul_x1: mul_in_matrix_pwas;
	signal signal_mul_x2: mul_in_matrix_pwas;
	signal signal_mul_y: mul_out_matrix_pwas;

signal start_sign : std_logic;

signal mulIn1_PWAS : mul_in_matrix_pwas;
signal mulIn2_PWAS : mul_in_matrix_pwas;
signal mulOut_PWAS : mul_out_matrix_pwas;

signal mulIn1_dyn : mul_in_matrix_pwas;
signal mulIn2_dyn : mul_in_matrix_pwas;
signal mulOut_dyn : mul_out_matrix_pwas;

signal pwas_start_signal : std_logic;
signal dynamic_signal : integer range 0 to N_DYN_PWAS-1;
signal valid_output_dyn : std_logic;
signal signal_stop : std_logic;
signal mux_sel : std_logic;
begin

	-- Port maps
	
	inst_pwasFunctionCircuit: discontinuousPwasFunctionCircuit
		port map ( clk => clk,
		         reset => reset,
		         start => pwas_start_signal,
					dynamic => dynamic_signal,
		             x => signal_x,
		         mul_y => mulOut_PWAS,
		        mul_x1 => mulIn1_PWAS,
		        mul_x2 => mulIn2_PWAS,
		             y => signal_y,
		          done => signal_stop);
					 
    findDynamic_portmap : findDynamic
    Port map( clk => clk,
           reset => reset,
			  start => start,
			  valid_output => valid_output_dyn,
			  mulIn1 => mulIn1_dyn,
			  mulIn2 => mulIn2_dyn,
			  mulOut => mulOut_dyn,
			  x => signal_x,
           dynamic => dynamic_signal);

	inst_mulBank : mulBank 
	port map (  x1 => signal_mul_x1,
		         x2 => signal_mul_x2,
		          y => signal_mul_y);

mul_process : process(signal_mul_y,mulIn1_PWAS,mulIn2_PWAS,mulIn1_dyn,mulIn2_dyn,mux_sel)
begin
if mux_sel = '1' then
        signal_mul_x1 <= mulIn1_PWAS;
        signal_mul_x2 <= mulIn2_PWAS;
else
        signal_mul_x1 <= mulIn1_dyn;
        signal_mul_x2 <= mulIn2_dyn;
end if;

    mulOut_PWAS <= signal_mul_y;
    mulOut_Dyn <= signal_mul_y;
end process;


mux_sel_process : process(clk,reset)
begin
if reset = '0' then
    mux_sel <= '0';
    pwas_start_signal <= '0';
elsif rising_edge(clk) then
    pwas_start_signal <= '0';
    if valid_output_dyn = '1' then
        mux_sel <= '1';
        pwas_start_signal <= '1';
    elsif signal_stop = '1' then
        mux_sel <= '0';
    end if;
end if;
end process;
			
	-- Assignments

done <= signal_stop;

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

end Behavioral;

