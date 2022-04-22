function writeOutputShift(object,ADCmeanValOut,cir_par)

% signed_ = true;
% if strcmp(cir_par.inputRepresentation,'unsigned')
%     signed_ = false;
% end

nx = object.getNumberOfStates;
nd = object.getNumberOfUnmeasurableInputs;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;

fout = fopen([cir_par.folder,'outputShift.vhdl'],'w');

fprintf(fout,'library IEEE;\n');
fprintf(fout,'use IEEE.STD_LOGIC_1164.all;\n');
fprintf(fout,'use IEEE.NUMERIC_STD.all;\n');
fprintf(fout,'use work.observer_package.all;\n\n');

fprintf(fout,'entity outputShift is\n');
fprintf(fout,'Port(');
fprintf(fout,'x_in : in compositeStateType_OBS;\n');
fprintf(fout,'	x_out : out compositeStateType_OBS');
fprintf(fout,');\n');
fprintf(fout,'end outputShift;\n\n');

fprintf(fout,'architecture Behavioral of outputShift is\n');
if nx+nd == 1
    fprintf(fout,'constant shiftConst : std_logic_vector(N_BIT_OBS downto 0) := ');
    num = decimal2signed(ADCmeanValOut(1),cir_par.inputResolution+1,0);
    fprintf(fout,'"%s";',num.bin);
else
    fprintf(fout,'type shiftConstType is array(0 to nx_OBS+nd_OBS-1) of std_logic_vector(N_BIT_OBS downto 0);\n');
    fprintf(fout,'constant shiftConst : shiftConstType := (');
    for i=1:nx+nd-1
        num = decimal2signed(ADCmeanValOut(i),cir_par.inputResolution+1,0);
        fprintf(fout,'"%s",',num.bin);
    end
    num = decimal2signed(ADCmeanValOut(end),cir_par.inputResolution+1,0);
    fprintf(fout,'"%s");\n',num.bin);
end

fprintf(fout,'begin\n');
fprintf(fout,'main_process : process(x_in)\n');
fprintf(fout,'variable tmp : std_logic_vector(N_BIT_OBS downto 0);\nbegin\n');

fprintf(fout,'for i in 0 to nx_OBS+nd_OBS-1 loop\n');
    fprintf(fout,'\ttmp := std_logic_vector(resize(signed(x_in(i)),N_BIT_OBS+1)+signed(shiftConst');
    if nx+nd ~= 1
        fprintf(fout,'(i)');
    end
    fprintf(fout,'));\n');
    fprintf(fout,'\tx_out(i) <= tmp(N_BIT_OBS-1 downto 0);\n');
fprintf(fout,'end loop;\n');

% fprintf(fout,'begin\n');
% fprintf(fout,'main_process : process(x_in)');
% fprintf(fout,'\nbegin\n');
% 
% fprintf(fout,'for i in 0 to nx_OBS+nd_OBS-1 loop\n');
% %     if signed_
% fprintf(fout,'\tx_out(i) <= std_logic_vector(resize(resize(signed(x_in(i)),N_BIT_OBS+1)+signed(shiftConst');
% if nx+nd ~= 1
%     fprintf(fout,'(i)');
% end
% fprintf(fout,'),N_BIT_OBS));\n');
% %     else
% %      fprintf(fout,'\tu_out(i) <= std_logic_vector(resize(signed(''0'' & u_in(i)),N_BIT_OBS+1)+signed(shiftConst(i)),N_BIT_OBS));\n');
% %     end
% fprintf(fout,'end loop;\n');

fprintf(fout,'end process;\n');
fprintf(fout,'end Behavioral;\n');