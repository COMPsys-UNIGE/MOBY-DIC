function writeInputShift(object,ADCmeanValIn,cir_par)

signed_ = true;
if strcmp(cir_par.inputRepresentation,'unsigned')
    signed_ = false;
end

nx = object.getNumberOfStates;
nd = object.getNumberOfUnmeasurableInputs;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;
nxref = numel(object.controller.getTrackingVariable);
fout = fopen([cir_par.folder,'inputShift.vhdl'],'w');

fprintf(fout,'library IEEE;\n');
fprintf(fout,'use IEEE.STD_LOGIC_1164.all;\n');
fprintf(fout,'use IEEE.NUMERIC_STD.all;\n');
fprintf(fout,'use work.embeddedSystemPackage.all;\n\n');

fprintf(fout,'entity inputShift is\n');
fprintf(fout,'Port(');
fprintf(fout,'u_in : in in_matrix_ES;\n');
fprintf(fout,'	u_out : out in_matrix_ES');
fprintf(fout,');\n');
fprintf(fout,'end inputShift;\n\n');

fprintf(fout,'architecture Behavioral of inputShift is\n');
if np+ny+nxref ~= 1
fprintf(fout,'type shiftConstType is array(0 to N_INPUT_ES-1) of std_logic_vector(N_BIT_ES downto 0);\n');
fprintf(fout,'constant shiftConst : shiftConstType := (');
for i=1:np+ny+nxref-1
    num = decimal2signed(-ADCmeanValIn(i),cir_par.inputResolution+1,0);
    fprintf(fout,'"%s",',num.bin);
end
num = decimal2signed(-ADCmeanValIn(end),cir_par.inputResolution+1,0);
fprintf(fout,'"%s");\n',num.bin);
else
    num = decimal2signed(-ADCmeanValIn,cir_par.inputResolution+1,0);
    fprintf(fout,'constant shiftConst : std_logic_vector(N_BIT_ES downto 0) := "%s";\n',num.bin);
end
fprintf(fout,'begin\n');
fprintf(fout,'main_process : process(u_in)');
fprintf(fout,'\nbegin\n');


    fprintf(fout,'for i in 0 to N_INPUT_ES-1 loop\n');
    if signed_
    fprintf(fout,'\tu_out(i) <= std_logic_vector(resize(resize(signed(u_in(i)),N_BIT_ES+1)+signed(shiftConst');
    if np+ny+nxref ~= 1
    fprintf(fout,'(i)');
    end
    fprintf(fout,'),N_BIT_ES));\n');
    else
     fprintf(fout,'\tu_out(i) <= std_logic_vector(resize(signed(''0'' & u_in(i))+signed(shiftConst');
    if np+ny+nxref ~= 1
    fprintf(fout,'(i)');
    end
    fprintf(fout,'),N_BIT_ES));\n');       
    end
fprintf(fout,'end loop;\n');

fprintf(fout,'end process;\n');
fprintf(fout,'end Behavioral;\n');