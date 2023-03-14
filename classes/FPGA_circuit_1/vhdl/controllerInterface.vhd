library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity controllerInterface is
Port ( clk : in STD_LOGIC;
	ce : in STD_LOGIC;
	reset : in STD_LOGIC;
	start : in STD_LOGIC;
	x1 : in STD_LOGIC_VECTOR (11 downto 0);
	x2 : in STD_LOGIC_VECTOR (11 downto 0);
	xref1 : in STD_LOGIC_VECTOR (11 downto 0);
	u1 : out STD_LOGIC_VECTOR (11 downto 0) := "011111111111";
	sample : out STD_LOGIC := '0';
	done : out STD_LOGIC := '0'
);
end controllerInterface;

architecture Behavioral of controllerInterface is

COMPONENT control
PORT (
	ap_clk : IN STD_LOGIC;
	ap_rst_n : IN STD_LOGIC;
	ap_start : IN STD_LOGIC;
	ap_done : OUT STD_LOGIC;
	ap_idle : OUT STD_LOGIC;
	ap_ready : OUT STD_LOGIC;
	x_in_0 : in STD_LOGIC_VECTOR (11 downto 0);
	x_in_1 : in STD_LOGIC_VECTOR (11 downto 0);
	ref_in : in STD_LOGIC_VECTOR (11 downto 0);
	u_opt : out STD_LOGIC_VECTOR (11 downto 0) := "011111111111";
	u_opt_ap_vld : OUT STD_LOGIC);
END COMPONENT;

signal admm_start : STD_LOGIC := '0';
signal admm_idle : STD_LOGIC;
signal admm_ready : STD_LOGIC;
signal u_opt_int : STD_LOGIC_VECTOR (11 downto 0) := "011111111111";
signal u_opt_vld : STD_LOGIC;
signal reset_n : STD_LOGIC;

begin

u1 <= "011111111111" when reset = '0' else 
	u_opt_int when u_opt_vld = '1' and rising_edge(clk);

admm_block : control
PORT MAP (
	ap_clk => clk,
	ap_rst_n => reset_n,
	ap_start => admm_start,
	ap_done => done,
	ap_idle => admm_idle,
	ap_ready => admm_ready,
	x_in_0 => x1,
	x_in_1 => x2,
	ref_in => xref1,
	u_opt => u_opt_int,
	u_opt_ap_vld => u_opt_vld);

process(clk, reset)
begin
	if reset = '0' then
		admm_start <= '0';
		sample <= '0';
	elsif rising_edge(clk) then
		if start = '1' then
			admm_start <= '1';
			sample <= '1';
		elsif admm_ready = '1' then
			admm_start <= '0';
			sample <= '0';
		else
			admm_start <= admm_start;
			sample <= '0';
		end if;
	end if;
end process;
end Behavioral;
