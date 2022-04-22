----------------------------------------------------------------------------------
-- MOBY-DIC Project
-- www.mobydic-project.eu
--
-- Project Name: pwag_par
--
-- Module Name: FSM
--
-- Description: 
-- This Finite State Machine corresponds to the binary search tree to be exlored 
-- for the computation of the PWAG function. Each state corresponds to a node of 
-- the tree
--
-- Copyright:
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

-- See file pwag_par_package.vhd for comments.
use work.pwagFunctionPackage.all;

-- Entity declaration

-- Inputs:
--  - clk: Clock signal
--  - reset: Reset signal (active low)
--  - decision: Indicates which branch of the tree must be followed

-- Outputs:
--  - leaf: If '1', the current state corresponds to a leaf node of the tree
--  - signal_addr: It is the memory address of the first coefficient related to the current state

entity FSM is
	Port ( clk : in std_logic;
		reset : in std_logic;
        restart : in std_logic;
	    decision : in std_logic;
		leaf : out std_logic;
		addr : out addr_array_PWAG;
        compute_done : in std_logic);
end FSM;

-- Architecture declaration

architecture behave of FSM is

	-- Signals declaration

	-- Indicates if a leaf node of the tree has been reached
	signal signal_leaf: std_logic;
	-- Signal defining the current and next state
	signal State, Next_State : unsigned(N_BIT_STATE_PWAG-1 downto 0);
	-- Addresses used to retrieve the coefficients from block Memory 
	signal signal_addr : addr_array_PWAG;

begin

	-- Next State combinational logic
	next_state_proc: process(State, decision, restart, signal_leaf)
	begin
		-- Until the computation of the MAC has not finished, or if a leaf node has been 
		-- reached, the current state is mantained
		if signal_leaf='1' then
            if restart = '1' then
                Next_State <= (others => '0');
            else
                Next_State<=State;
            end if;
		else
			case State is
			
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---
				
				when others => Next_State <= (others => '0');
			end case;
		end if;
	end process;

	--State register
	reg_proc : process(clk,reset)
	begin
		if(reset = '0') then
			State <= (others=>'0');
		elsif rising_edge(clk) then
            if compute_done = '1'  or signal_leaf = '1'  then
                State <= Next_State;
            end if;
        end if;
	end process;

	-- Output combinational logic
	-- This process gives the addresses associated to the coefficients of the different
	-- output functions. This obviously makes sense only for leaf nodes: for non-leaf nodes 
	-- only signal_addr(0) is used, therefore all other signals are assigned as zeros in order 
	-- to improve synthesis. 
	output_proc : process(State)
	begin
	--Set output as function of current state and imput
		case State is
		
--- BEGIN MATLABGEN ---
--- END MATLABGEN ---
						
		end case;
	end process;

	-- If reset is 1, the matrix of addresses is put in output
	addr_proc: process(signal_addr,reset)
		begin
			if reset='0' then 
				for i in 0 to N_FUN_PWAG-1 loop
					addr(i) <= (others => '0');
				end loop;
			else
				addr <= signal_addr;
			end if;
		end process;
		

	-- If reset is 0, output leaf is 0
	leaf <= signal_leaf and reset;
	
end behave;
