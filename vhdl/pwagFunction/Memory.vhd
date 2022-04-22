----------------------------------------------------------------------------------
-- MOBY-DIC Project
-- www.mobydic-project.eu
--
-- Project Name: pwag_par
--
-- Module Name: Memory
--
-- Description: 
-- Block containing as many LUTs as the number of the domain dimensions plus one.
-- LUTs L1, ..., Ln and R contain the coefficients of the affine expression
-- L1 x1 + L2 x2 + ... + Ln xn + R
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

-- See file pwag_ser_package.vhd for comments.
use work.pwagFunctionPackage.all;

-- Entity declaration

-- Inputs:
--  - addr: Address of the memory entry
--  - clk: Clock signal

-- Outputs:
--  - res: Value of the memory entries selected with addr

entity Memory is
	Port ( addr : in unsigned(N_BIT_ADDR_PWAG-1 downto 0);
		clk : in std_logic;
		res : out memory_out_matrix_PWAG);
end Memory;

-- Architecture declaration

architecture Behavioral of Memory is

--- BEGIN MATLABGEN ---
--- END MATLABGEN ---

begin

read_process : process(clk)
begin
if rising_edge(clk) then
    for i in 0 to N_DIM_PWAG loop
        res(i) <= ROM(to_integer(addr))(i);
    end loop;
end if;
end process;
	
end Behavioral;
