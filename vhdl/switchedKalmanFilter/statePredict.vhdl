----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    16:51:51 01/04/2016 
-- Design Name: 
-- Module Name:    nextStateBlock - Behavioral 
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

entity statePredict is
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
end statePredict;

architecture Behavioral of statePredict is

component shifterMACPredictor is Port ( reset : in  STD_LOGIC;
           clk : in  STD_LOGIC;
			  clr : in std_logic;
			  mulOut : in std_logic_vector(N_BIT_COEFF_OBS+N_BIT_OBS-1 downto 0);
			  nShift : in integer range 0 to nBusMacPredict_OBS-1;
           output : out  std_logic_vector(nBusMacPredict_OBS-1 downto 0));
end component;

type input_merge_type is array(nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);
type outputMac_array is array(nx_OBS+nd_OBS-1 downto 0) of std_logic_vector(nBusMacPredict_OBS-1 downto 0);

signal internal_counter : integer range 0 to nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS+4;
signal input_merge : input_merge_type;
signal MAC_clr : std_logic;
signal index_counter : integer range 0 to nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS;
signal shiftI : integer range 0 to nBusMacPredict_OBS-1;
signal macOut : outputMac_array;
signal shift_counter : integer range 0 to nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS;
begin

-- Internal counter process
internal_counter_process : process(reset,clk)
begin
if reset = '0' then
	internal_counter <= 0;
elsif rising_edge(clk) then
	if internal_counter = nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS+4 then 
		internal_counter <= 0;
	else
		if internal_counter /= 0 or start = '1' then 
			internal_counter <= internal_counter +1;
		end if;
	end if;
end if;
end process;

-- index_counter is equal to internal counter saturate to nx+nd+np+nu+ny
index_counter <= internal_counter-1 when (internal_counter < nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS+2 and internal_counter > 0)
						else 0;

shift_counter <= internal_counter-2 when (internal_counter < nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS+3 and internal_counter > 1)
						else 0;

--signal to multiplier process
sign_to_mul_process : process(input_merge,index_counter,dynamic)
begin
for i in 0 to nx_OBS+nd_OBS-1 loop
	mulIn1(i) <= resize(signed(input_merge(index_counter)),N_BIT_MUL_OBS);
	mulIn2(i) <= resize(signed(predictMatrix(dynamic)(i)(index_counter)),N_BIT_MUL_OBS);
end loop;
end process;

shiftI <= preShiftPredict_OBS(dynamic)(shift_counter);

-- MAC port map
MAC_GENERATE: for I in 0 to nx_OBS+nd_OBS-1 generate
      MACX : shifterMACPredictor port map (reset => reset,
			clk => clk,
			clr => MAC_clr,
			mulOut => std_logic_vector(mulOut(I)(N_BIT_OBS+N_BIT_COEFF_OBS-1 downto 0)),
			nShift => shiftI,
			output => macOut(I));
   end generate MAC_GENERATE;

-- MAC out process
macOut_process: process(macOut,dynamic)
variable tmp : std_logic_vector(N_BIT_COEFF_OBS+nBusMacPredict_OBS-1 downto 0); 
begin
	for i in 0 to nx_OBS+nd_OBS-1 loop
		tmp := std_logic_vector(signed(alphaPred_OBS(dynamic))*signed(macOut(i))+signed(betaPred_OBS(dynamic)));
		x_next(i) <= tmp(nDecOutPred_OBS(dynamic)+N_BIT_OBS-1 downto nDecOutPred_OBS(dynamic));
	end loop;
end process;

-- MAC_clr assign
MAC_clr_process : process(internal_counter,start)
begin
MAC_Clr <= '0';
if internal_counter = 0 then
	MAC_Clr <= '1';
end if;
end process;

-- valid output assign
with internal_counter select
	valid_output <= '1' when  nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS+4,
					    '0' when others;

-- assign input
input_process : process(x,u)
begin
for i in 0 to nx_OBS+nd_OBS-1 loop
	input_merge(i) <= x(i);
end loop;
for i in 0 to np_OBS+nu_OBS+ny_OBS-1 loop
	input_merge(i+nx_OBS+nd_OBS) <= u(i);
end loop;
input_merge(nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS) <= inputOne_OBS;
end process;



end Behavioral;

