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

entity pwagFunctionCircuit is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	     start : in std_logic;
	         x : in x_matrix_PWAG;
	     mul_y : in mul_out_matrix_PWAG;
	    mul_x1 : out mul_in_matrix_PWAG;
	    mul_x2 : out mul_in_matrix_PWAG;
	         y : out y_matrix_PWAG;
	      done : out  std_logic);
end pwagFunctionCircuit;

-- Architecture declaration

architecture Behavioral of pwagFunctionCircuit is

component compute is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	      leaf : in std_logic;
	   restart : in std_logic;
		   addr : in addr_array_PWAG;
	     mul_y : in mul_out_matrix_PWAG;
	    mul_x2 : out mul_in_matrix_PWAG;
	         y : out y_matrix_PWAG;
      decision: out std_logic;
	      compute_done: out  std_logic;
			compute_func : out std_logic);
end component;

component FSM is
	Port ( clk : in std_logic;
		reset : in std_logic;
        restart : in std_logic;
	    decision : in std_logic;
		leaf : out std_logic;
		addr : out addr_array_PWAG;
        compute_done : in std_logic);
end component;


signal addr_from_FSM : addr_array_PWAG;

signal signal_decision : std_logic;
signal signal_leaf : std_logic;

signal compute_done : std_logic;
signal function_done : std_logic;
begin	

compute_portmap : compute 
port map ( clk => clk,
	     reset => reset,
	     leaf => signal_leaf,
		  addr => addr_from_FSM,
		  restart => start,
	     mul_y => mul_y,
	    mul_x2 => mul_x2,
	         y => y,
      decision => signal_decision,
	      compute_done => compute_done,
			compute_func => function_done);


FSM_portmap : FSM
port map ( clk => clk,
		reset => reset,
        restart => start,
	    decision => signal_decision,
		leaf => signal_leaf,
		addr => addr_from_FSM,
		compute_done => compute_done);



mul_signal_process : process(x)
begin
for i in 0 to N_DIM_PWAG-1 loop
	mul_x1(i) <= resize(x(i),N_BIT_MUL_PWAG);
end loop;
end process;

done <= function_done;

end Behavioral;