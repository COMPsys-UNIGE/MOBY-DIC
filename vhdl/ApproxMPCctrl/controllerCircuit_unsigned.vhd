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
use work.pwasFunctionPackage.all;
-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity controllerCircuit is
	Port( clk : in std_logic;
			reset : in std_logic;
			start : in std_logic;
			 done : out std_logic;
                x : in x_matrix_CTRL;
            mul_y : in mul_out_matrix_CTRL;
           mul_x1 : out mul_in_matrix_CTRL;
           mul_x2 : out mul_in_matrix_CTRL;
                u : out u_matrix_CTRL
);
end controllerCircuit;

architecture Behavioral of controllerCircuit is

component pwasFunctionCircuit is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	     start : in std_logic;
	         x : in x_matrix_PWAS;
	     mul_y : in mul_out_matrix_PWAS;
	    mul_x1 : out mul_in_matrix_PWAS;
	    mul_x2 : out mul_in_matrix_PWAS;
	         y : out y_matrix_PWAS;
	      done : out  std_logic);
end component;

-- Signals declaration
	
	signal signal_y: y_matrix_PWAS;
	signal signal_x: x_matrix_PWAS;
	signal signal_mul_y : mul_out_matrix_PWAS;
	signal signal_mul_x1 :  mul_in_matrix_PWAS;
	signal signal_mul_x2 :  mul_in_matrix_PWAS;


begin

	-- Port maps
	
	inst_pwasFunctionCircuit: pwasFunctionCircuit
		port map ( clk => clk,
		         reset => reset,
		         start => start,
		             x => signal_x,
		         mul_y => signal_mul_y,
		        mul_x1 => signal_mul_x1,
		        mul_x2 => signal_mul_x2,
		             y => signal_y,
		          done => done);

-- signal assign
main_proc : process(x,mul_y,signal_mul_x1,signal_mul_x2,signal_y)
begin
for i in 0 to N_DIM_PWAS loop
     mul_x1(i) <= signal_mul_x1(i);
     mul_x2(i) <= signal_mul_x2(i);

	  signal_mul_y(i) <= mul_y(i);
end loop;
for i in 0 to N_DIM_PWAS-1 loop
     signal_x(i) <= "0" & x(i);
end loop;
for i in 0 to N_FUN_PWAS-1 loop
     u(i) <= signal_y(i);
end loop;
end process;


end Behavioral;

