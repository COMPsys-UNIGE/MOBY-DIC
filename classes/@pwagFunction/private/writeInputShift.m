function writeInputShift(object,ADCmeanValIn,cir_par)

signed_ = true;
if strcmp(cir_par.inputRepresentation,'unsigned')
    signed_ = false;
end

ndim = object.getDomainDimensions;
nx = numel(ADCmeanValIn);
fout = fopen([cir_par.folder,'inputShift.vhdl'],'w');

fprintf(fout,'library IEEE;\n');
fprintf(fout,'use IEEE.STD_LOGIC_1164.all;\n');
fprintf(fout,'use IEEE.NUMERIC_STD.all;\n');
fprintf(fout,'use work.pwagFunctionPackage.all;\n\n');

fprintf(fout,'entity inputShift is\n');
fprintf(fout,'Port(');

fprintf(fout,'x_in : in x_matrix_PWAG;\n');
fprintf(fout,'	x_out : out x_matrix_PWAG);\n');

fprintf(fout,'end inputShift;\n\n');

fprintf(fout,'architecture Behavioral of inputShift is\n');
if ndim ~= 1
    fprintf(fout,'type shiftConstType is array(0 to N_DIM_PWAG-1) of signed(N_BIT_PWAG downto 0);\n');
    fprintf(fout,'constant shiftConst : shiftConstType := (');
    for i=1:nx-1
        num = decimal2signed(-ADCmeanValIn(i),cir_par.inputResolution+1,0);
        fprintf(fout,'"%s",',num.bin);
    end
    num = decimal2signed(-ADCmeanValIn(end),cir_par.inputResolution+1,0);
    fprintf(fout,'"%s");\n',num.bin);
else
    numi = decimal2signed(-ADCmeanValIn,cir_par.inputResolution+1,0);
    fprintf(fout,'constant shiftConst : signed(N_BIT_PWAG downto 0) := "%s";\n',numi.bin);
end
fprintf(fout,'begin\n');
fprintf(fout,'main_process : process(x_in)\nbegin\n');

fprintf(fout,'for i in 0 to N_DIM_PWAG-1 loop\n');
if signed_
    fprintf(fout,'\tx_out(i) <= resize(resize(signed(x_in(i)),N_BIT_PWAG+1)+shiftConst');
    if ndim ~= 1
        fprintf(fout,'(i)');
    end
    fprintf(fout,',N_BIT_PWAG);\n');
else
    fprintf(fout,'\tx_out(i) <= resize(signed(''0'' & x_in(i))+shiftConst');
    if ndim ~= 1
        fprintf(fout,'(i)');
    end
    fprintf(fout,',N_BIT_PWAG);\n');
end
fprintf(fout,'end loop;\n');

fprintf(fout,'end process;\n');
fprintf(fout,'end Behavioral;\n');