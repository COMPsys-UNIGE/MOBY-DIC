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

entity pwasFunctionCircuit is
	Port ( clk : in std_logic;
	     reset : in std_logic;
	     start : in std_logic;
	         x : in x_matrix_pwas;
	     mul_y : in mul_out_matrix_pwas;
	    mul_x1 : out mul_in_matrix_pwas;
	    mul_x2 : out mul_in_matrix_pwas;
	         y : out y_matrix_pwas;
	      done : out  std_logic);
end pwasFunctionCircuit;

-- Architecture declaration

architecture Behavioral of pwasFunctionCircuit is

	-- Signals declaration

	signal signal_mul_y : mul_out_matrix_pwas;
	signal signal_scaleA : coeff_matrix_pwas;
	signal signal_x : x_matrix_pwas;
	signal signal_x_int : int_matrix_pwas;
	signal signal_x_dec : dec_matrix_pwas;
	signal signal_x_dec_sorted : dec_matrix_pwas;
	signal signal_mu : mu_matrix_pwas;
	signal signal_address : addr_matrix_pwas;
	signal signal_w: mem_matrix_pwas;
	signal signal_mul_done : std_logic;
	signal signal_mul_done_scale : std_logic;
	signal signal_sorter_done : std_logic;
	signal signal_sorter_start : std_logic;
	signal signal_mu_done : std_logic;
	signal signal_address_done : std_logic;
	signal signal_mem_done : std_logic;
	signal signal_mul_done_final : std_logic;
	signal signal_mul_start : std_logic;
	signal signal_mu_start : std_logic;
	signal signal_address_start : std_logic;
	signal signal_mem_start : std_logic;
	signal signal_sum_start : std_logic;
	signal signal_scale_start : std_logic;
	signal signal_scale_done : std_logic;
		
	-- Components declaration
	 
	component mulBankManager is
		Port (     clk : in std_logic;
		         reset : in std_logic;
		         start : in std_logic;
		   sorter_done : in std_logic;
		      mem_done : in std_logic;
		             x : in x_matrix_pwas;
		        scaleA : in coeff_matrix_pwas;
		            mu : in mu_matrix_pwas;
		             w : in mem_matrix_pwas;
		      mul_done : in std_logic;
		mul_done_scale : out std_logic;
		mul_done_final : out std_logic;
		        mul_x1 : out mul_in_matrix_pwas;
		        mul_x2 : out mul_in_matrix_pwas;
		     mul_start : out std_logic);
	end component;
	
	component Scale is
		Port ( clk : in std_logic;
	        reset : in std_logic;
	        start : in std_logic;
	            x : in x_matrix_pwas;
	         Axmb : in mul_out_matrix_pwas;
	          xmb : out x_matrix_pwas;
	       scaleA : out coeff_matrix_pwas;
	        x_int : out int_matrix_pwas;
	        x_dec : out dec_matrix_pwas;
	         done : out  std_logic);
	end component;

	component Sorter is
		Port ( clk : in std_logic;
		     start : in std_logic;
		     reset : in std_logic;
		         x : in dec_matrix_pwas;
		         y : out dec_matrix_pwas;
		      done : out std_logic);
	end component;
	
	component mu_generator is
		Port (   clk : in std_logic;
		       reset : in std_logic;
		       start : in std_logic;
		x_dec_sorted : in dec_matrix_pwas;
		          mu : out mu_matrix_pwas;
		        done : out std_logic);
	end component;
	
	component address_generator is
		Port (   clk : in std_logic;
		       reset : in std_logic;
		       start : in std_logic;
		       x_int : in int_matrix_pwas;
		       x_dec : in dec_matrix_pwas;
		x_dec_sorted : in dec_matrix_pwas;
		     address : out addr_matrix_pwas;
		        done : out std_logic);
	end component;
	
	component Memory is
		Port ( clk : in std_logic;
		     start : in std_logic;
		     reset : in std_logic;
		      addr : in addr_matrix_pwas;
		       res : out mem_matrix_pwas;
		      done : out std_logic);
	end component;
	
	component Sum is
		Port ( clk : in  std_logic;
		     reset : in  std_logic;
		     start : in  std_logic;
		         x : in  mul_out_matrix_pwas;
		         y : out y_matrix_pwas;
		      done : out  std_logic);
	end component;
			
begin

	inst_mulBankManager : mulBankManager
	port map ( clk => clk,
	         reset => reset,
	         start => start,
	   sorter_done => signal_sorter_done,
	      mem_done => signal_mem_done,
	             x => signal_x,
	        scaleA => signal_scaleA,
	            mu => signal_mu,
	             w => signal_w,
	      mul_done => signal_mul_done,
	mul_done_scale => signal_mul_done_scale,
	mul_done_final => signal_mul_done_final,
	        mul_x1 => mul_x1,
	        mul_x2 => mul_x2,
	     mul_start => signal_mul_start);
		  
	inst_scale: Scale
	port map ( clk => clk,
	         reset => reset,
	         start => signal_scale_start,
	             x => x,
	          Axmb => signal_mul_y,
				  xmb => signal_x,
	        scaleA => signal_scaleA,
	         x_int => signal_x_int,
	         x_dec => signal_x_dec,
	          done => signal_scale_done);	 

	 inst_sorter : Sorter
	port map ( clk => clk,
	         reset => reset,
	         start => signal_sorter_start,
	             x => signal_x_dec,
	             y => signal_x_dec_sorted,
	          done => signal_sorter_done);
	 
	inst_mu_generator : mu_generator
	port map ( clk => clk,
	         reset => reset,
	         start => signal_mu_start,
	  x_dec_sorted => signal_x_dec_sorted,
	            mu => signal_mu,
	          done => signal_mu_done);	 

	inst_address_generator : address_generator
	port map ( clk => clk,
	         reset => reset,
	         start => signal_address_start,
	         x_int => signal_x_int,
	         x_dec => signal_x_dec,
	  x_dec_sorted => signal_x_dec_sorted,
	       address => signal_address,
	          done => signal_address_done);
	 
   inst_Memory : Memory
	port map ( clk => clk,
	         reset => reset,
	         start => signal_mem_start,
	          addr => signal_address,
	           res => signal_w,
	          done => signal_mem_done); 
	 
	inst_Sum : Sum
	port map ( clk => clk,
	         reset => reset,
	         start => signal_sum_start,
	             x => signal_mul_y,
	             y => y,
	          done => done);
	 
	signal_sorter_start <= signal_scale_done;
	signal_mu_start <= signal_sorter_done;	
	signal_address_start <= signal_sorter_done;	
	signal_mem_start <= signal_address_done;	
	signal_sum_start <= signal_mul_done_final;					 
	signal_scale_start <= signal_mul_done_scale;
	
	proc_mul : process(clk,reset,signal_mul_start)
	begin
		if reset = '0' then
			for i in 0 to N_DIM_PWAS loop
				signal_mul_y(i) <= (others => '0');
			end loop;
		elsif rising_edge(clk) and signal_mul_start = '1' then
			signal_mul_y <= mul_y;
		end if;
	end process;
		
	signal_mul_done <= '0' when reset = '0' else
	                   signal_mul_start when rising_edge(clk);
	
end Behavioral;

