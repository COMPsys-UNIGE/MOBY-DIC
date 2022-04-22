----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    08:59:44 01/15/2016 
-- Design Name: 
-- Module Name:    mulBank - Behavioral 
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

use work.pwagFunctionPackage.all;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity mulBank is
Port(			  x1 : in mul_in_matrix_PWAG;
			  x2 : in mul_in_matrix_PWAG;
			  y : out mul_out_matrix_PWAG);
end mulBank;

architecture Behavioral of mulBank is

begin

main_mul_process : process(x1,x2)
begin
for i in 0 to N_DIM_PWAG-1 loop
	y(i) <= x1(i)*x2(i);
end loop;
end process;

end Behavioral;

