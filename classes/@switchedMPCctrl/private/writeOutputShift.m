function writeOutputShift(object,ADCmeanValOut,cir_par)

% signed_ = true;
% if strcmp(cir_par.inputRepresentation,'unsigned')
%     signed_ = false;
% end

nout = numel(ADCmeanValOut);

fout = fopen([cir_par.folder,'outputShift.vhdl'],'w');

fprintf(fout,'library IEEE;\n');
fprintf(fout,'use IEEE.STD_LOGIC_1164.all;\n');
fprintf(fout,'use IEEE.NUMERIC_STD.all;\n');
fprintf(fout,'use work.controller_package.all;\n\n');

fprintf(fout,'entity outputShift is\n');
fprintf(fout,'Port(');

    fprintf(fout,'u_in : in u_matrix_CTRL;\n');
    fprintf(fout,'	u_out : out u_matrix_CTRL');
fprintf(fout,');\n');
fprintf(fout,'end outputShift;\n\n');

fprintf(fout,'architecture Behavioral of outputShift is\n');
if nout ~= 1
fprintf(fout,'type shiftConstType is array(0 to N_FUN_CTRL-1) of signed(N_BIT_OUT_CTRL downto 0);\n');
fprintf(fout,'constant shiftConst : shiftConstType := (');
for i=1:nout-1
    num = decimal2signed(ADCmeanValOut(i),cir_par.outputResolution+1,0);
    fprintf(fout,'"%s",',num.bin);
end
num = decimal2signed(ADCmeanValOut(end),cir_par.outputResolution+1,0);
fprintf(fout,'"%s");\n',num.bin);
else
    numi = decimal2signed(ADCmeanValOut,cir_par.outputResolution+1,0);
    fprintf(fout,'constant shiftConst : signed(N_BIT_OUT_CTRL downto 0) := "%s";\n',numi.bin);
end
fprintf(fout,'begin\n');
fprintf(fout,'main_process : process(u_in)\n');
fprintf(fout,'variable tmp : std_logic_vector(N_BIT_OUT_CTRL downto 0);\nbegin\n');

fprintf(fout,'for i in 0 to N_FUN_CTRL-1 loop\n');
    fprintf(fout,'\ttmp := std_logic_vector(resize(signed(u_in(i)),N_BIT_OUT_CTRL+1)+shiftConst');
    if nout ~= 1
        fprintf(fout,'(i)');
    end
    fprintf(fout,');\n');
    fprintf(fout,'\tu_out(i) <= signed(tmp(N_BIT_OUT_CTRL-1 downto 0));\n');
fprintf(fout,'end loop;\n');



fprintf(fout,'end process;\n');
fprintf(fout,'end Behavioral;\n');