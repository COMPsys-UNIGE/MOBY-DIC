----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    17:10:56 01/04/2016 
-- Design Name: 
-- Module Name:    MAC - Behavioral 
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

use work.observer_package.all;
use IEEE.STD_LOGIC_1164.ALL;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity shifterMACPredictor is
    Port ( reset : in  STD_LOGIC;
           clk : in  STD_LOGIC;
			  clr : in std_logic;
			  mulOut : in std_logic_vector(N_BIT_COEFF_OBS+N_BIT_OBS-1 downto 0);
			  nShift : in integer range 0 to nBusMacPredict_OBS-1;
           output : out  std_logic_vector(nBusMacPredict_OBS-1 downto 0));
end shifterMACPredictor;

architecture Behavioral of shifterMACPredictor is
signal acc : std_logic_vector(nBusMacPredict_OBS-1 downto 0);
signal next_acc : std_logic_vector(nBusMacPredict_OBS-1 downto 0);
signal inreg : std_logic_vector(N_BIT_COEFF_OBS+N_BIT_OBS-1 downto 0);
signal reg_interm : std_logic_vector(nBusMacPredict_OBS-1 downto 0);
begin

in_reg_process : process(reset,clk)
begin
if reset = '0' then
	inreg <= (others => '0');
elsif rising_edge(clk) then
	if clr = '0' then
	inreg <= mulOut;
	else
		inreg <= (others => '0');
	end if;
end if;
end process;

inter_reg_process : process(reset,clk)
begin
if reset = '0' then
	reg_interm <= (others => '0');
elsif rising_edge(clk) then
	if clr = '0' then
		reg_interm <= std_logic_vector(SHIFT_LEFT(resize(signed(inreg),nBusMacPredict_OBS),nShift) );
	else
		reg_interm <= (others => '0');		
	end if;
end if;
end process;

reg_proc : process(reset,clk)
begin
if reset = '0' then
	acc <= (others => '0');
elsif rising_edge(clk) then
	acc <= next_acc;
end if;
end process;

--tmp <= signed(in1)*signed(in2);

next_acc <= std_logic_vector(signed(reg_interm) + signed(acc)) when clr = '0'
								else (others => '0');

output <= acc;

end Behavioral;

