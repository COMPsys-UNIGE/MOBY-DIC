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
use work.pwasFunctionPackage.all;

entity bitonic4 is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	     start : in std_logic;
	        x1 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        x2 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        x3 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        x4 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        up : in std_logic;
	        y1 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        y2 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        y3 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	        y4 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	      done : out std_logic);
end bitonic4;

architecture Behavioral of bitonic4 is

	signal state : integer range 0 to 3 := 0;
	signal next_state : integer range 0 to 3 := 0;
	signal comp1_x1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp1_x2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp1_y1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp1_y2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp1_up : std_logic;
	signal comp1_start : std_logic;
	signal comp2_x1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp2_x2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp2_y1 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp2_y2 : unsigned(N_BIT_COEFF_PWAS-1 downto 0);
	signal comp2_up : std_logic;
	signal comp2_start : std_logic;
	
	component comparator is
		Port ( clk : in std_logic;
		     start : in std_logic;
		     reset : in std_logic;
		        x1 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        x2 : in unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        y1 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        y2 : out unsigned(N_BIT_COEFF_PWAS-1 downto 0);
		        up : in std_logic);
	end component;

begin

	comp1: comparator
		port map( clk => clk,
		        reset => reset,
		        start => comp1_start,
		           x1 => comp1_x1,
		           x2 => comp1_x2,
		           y1 => comp1_y1,
		           y2 => comp1_y2,
		           up => comp1_up);

	comp2: comparator
		port map( clk => clk,
		        reset => reset,
		        start => comp2_start,
		           x1 => comp2_x1,
		           x2 => comp2_x2,
		           y1 => comp2_y1,
		           y2 => comp2_y2,
		           up => comp2_up);


-- Next State combinational logic
	next_state_proc: process(state,start)
	begin
		if (state = 0) then
			if start <= '0' then
				next_state <= 0;
			else
				next_state <= 1;
			end if;
		elsif (state = 3) then
			next_state <= 0;
		else
			next_state <= state+1;
		end if;
	end process;

	--State register
	reg_proc : process(clk,reset)
	begin
		if(reset = '0') then
			state <= 0;
		elsif rising_edge(clk) then
			state <= next_state;
		end if;
	end process;

	-- Output combinational logic
	output_proc : process(State,x1,x2,x3,x4,comp1_y1,comp1_y2,comp2_y1,comp2_y2,up,start)
	begin
		--Set output as function of current state and imput
		case state is
			when 0 =>
				comp1_x1 <= x1;
				comp1_x2 <= x2;
				comp1_up <= '0';
				comp1_start <= start;
				comp2_x1 <= x3;
				comp2_x2 <= x4;
				comp2_up <= '1';
				comp2_start <= start;
				done <= '0';
			when 1 =>
				comp1_x1 <= comp1_y1;
				comp1_x2 <= comp2_y1;
				comp1_up <= up;
				comp1_start <= '1';
				comp2_x1 <= comp1_y2;
				comp2_x2 <= comp2_y2;
				comp2_up <= up;
				comp2_start <= '1';
				done <= '0';
			when 2 =>
				comp1_x1 <= comp1_y1;
				comp1_x2 <= comp2_y1;
				comp1_up <= up;
				comp1_start <= '1';
				comp2_x1 <= comp1_y2;
				comp2_x2 <= comp2_y2;
				comp2_up <= up;
				comp2_start <= '1';
				done <= '0';
			when 3 =>
				comp1_x1 <= comp1_y1;
				comp1_x2 <= comp2_y1;
				comp1_up <= up;
				comp1_start <= '0';
				comp2_x1 <= comp1_y2;
				comp2_x2 <= comp2_y2;
				comp2_up <= up;
				comp2_start <= '0';
				done <= '1';
		end case;
	end process;

	y1 <= comp1_y1;
	y2 <= comp1_y2;
	y3 <= comp2_y1;
	y4 <= comp2_y2;


end Behavioral;

