function generatePackage(object,folder,FintPred, GintPred, alphaPred, betaPred, nbit, nint, nbit_coeff, nint_coeffPred,initialState,sampling_latency)

ndec = nbit-nint;
ndec_coeffPred = nbit_coeff-nint_coeffPred;


shiftMulPred = [ndec 0]+ndec_coeffPred;
preshiftPred = max(shiftMulPred)-shiftMulPred;

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;

M = [FintPred GintPred];

fout = fopen([folder,'observer_package.vhdl'],'w');

fprintf(fout,'library IEEE;\n');
fprintf(fout,'use IEEE.STD_LOGIC_1164.all;\n\n');
fprintf(fout,'use IEEE.NUMERIC_STD.all;\n\n');
fprintf(fout,'package observer_package is\n\n');

fprintf(fout,'constant N_BIT_OBS : integer := %d;\n',nbit);
fprintf(fout,'constant N_BIT_COEFF_OBS : integer := %d;\n\n',nbit_coeff);

fprintf(fout,'constant sampling_latency_OBS : integer := %d;\n\n',fix(sampling_latency));

fprintf(fout,'constant nx_OBS : integer := %d;\n',object.getNumberOfStates);
fprintf(fout,'constant np_OBS : integer := %d;\n',object.getNumberOfParameters);
fprintf(fout,'constant nd_OBS : integer := %d;\n',object.getNumberOfUnmeasurableInputs);
fprintf(fout,'constant nu_OBS : integer := %d;\n',object.getNumberOfInputs);
fprintf(fout,'constant ny_OBS : integer := %d;\n\n',object.getNumberOfOutputs);

% if nx ~= 0
%     fprintf(fout,'type xType_OBS is array(nx_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n');
% end
% if nd ~= 0
%     fprintf(fout,'type dType_OBS is array(nd_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n');
% end
% if nu ~= 0
%     fprintf(fout,'type uType_OBS is array(nu_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n');
% end
% if np ~= 0
%     fprintf(fout,'type pType_OBS is array(np_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n');
% end
% if ny ~= 0
%     fprintf(fout,'type yType_OBS is array(ny_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n\n');
% end

fprintf(fout,'constant minShiftPred_OBS : integer := %d;\n\n',min(preshiftPred));
fprintf(fout,'constant maxShiftPred_OBS : integer := %d;\n\n',max(preshiftPred));
uno = decimal2signed(1,nbit,0);
nBusMacPredict = nbit+nbit_coeff+ceil(log2(nx+np+nd+nu+ny))+max(preshiftPred)-min(preshiftPred);
fprintf(fout,'constant nBusMacPredict_OBS : integer := %d;\n',nBusMacPredict);%N_BIT_OBS+N_BIT_COEFF_OBS+nx_OBS+np_OBS+nd_OBS+nu_OBS+ny_OBS+maxShiftPred_OBS-minShiftPred_OBS;\n\n');
fprintf(fout,'constant inputOne_OBS : std_logic_vector(N_BIT_OBS-1 downto 0) := "%s"; -- uno messo bene\n\n',uno.bin);

fprintf(fout,'type compositeStateType_OBS is array (nx_OBS+nd_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n');
fprintf(fout,'type compositeInputType_OBS is array (np_OBS+nu_OBS+ny_OBS-1 downto 0) of std_logic_vector(N_BIT_OBS-1 downto 0);\n\n');
fprintf(fout,'\tconstant initialState : compositeStateType_OBS := (');
for i = 1:nx+nd-1
    numi = decimal2signed(initialState(nx+nd+1-i),nbit,ndec(nx+nd+1-i));
    fprintf(fout,'"%s",',numi.bin);
end
numi = decimal2signed(initialState(1),nbit,ndec(1));
fprintf(fout,'"%s");\n',numi.bin);

fprintf(fout,'\n');
fprintf(fout,'type observerMatrixColumn_OBS is array (0 to nx_OBS+nd_OBS+np_OBS+nu_OBS+ny_OBS) of std_logic_vector(N_BIT_COEFF_OBS-1 downto 0);\n');
fprintf(fout,'type observerMatrixType_OBS is array (0 to nx_OBS+nd_OBS-1) of observerMatrixColumn_OBS;\n\n');

alphaPredBin = decimal2signed(alphaPred,nbit_coeff);
betaPredBin = decimal2signed(betaPred,nBusMacPredict+nbit_coeff,max(shiftMulPred)+alphaPredBin.pointposition);


fprintf(fout,'constant alphaPred_OBS : std_logic_vector(N_BIT_COEFF_OBS-1 downto 0) := "%s";\n\n',alphaPredBin.bin);
%fprintf(fout,'constant betaPred : std_logic_vector(N_BIT_COEFF+N_BIT-1 downto 0) := "%s";\n\n',betaPredBin.bin);
fprintf(fout,'constant betaPred_OBS : std_logic_vector(N_BIT_COEFF_OBS+nBusMacPredict_OBS-1 downto 0) := "%s";\n\n',betaPredBin.bin);

% Write M matrix
fprintf(fout,'constant predictMatrix : observerMatrixType_OBS := (( ');
for i=1:size(M,1)
    for j=1:size(M,2)
        num = decimal2signed(M(i,j),nbit_coeff,0);
        fprintf(fout,'"%s"',num.bin);
        if j ~= size(M,2)
            fprintf(fout,',');
        else
            fprintf(fout,')');
        end
    end
    if i ~= size(M,1)
        fprintf(fout,',\n(');
    else
        fprintf(fout,');\n\n');
    end
end

fprintf(fout,'\ttype int_shift_array_OBS is array(0 to nx_OBS+nd_OBS+nu_OBS+np_OBS+ny_OBS) of integer;\n');
fprintf(fout,'\ttype int_array_OBS is array(0 to nx_OBS+nd_OBS+nu_OBS+np_OBS+ny_OBS-1) of integer;\n');

fprintf(fout,'\tconstant nDecOutPred_OBS : integer := %d;\n',alphaPredBin.pointposition+max(shiftMulPred));
%fprintf(fout,'\tconstant nDecMacOutPred_OBS : integer := %d;\n',max(shiftMulPred));
% fprintf(fout,'\tconstant nDec_OBS : int_array_OBS := (');

% for i = 1:nx+nd+np+nu+ny-1
%     fprintf(fout,'%d,',ndec(i));
% end
% fprintf(fout,'%d);\n',ndec(nx+nd+np+nu+ny));

fprintf(fout,'\n');
fprintf(fout,'\tconstant preShiftPredict_OBS : int_shift_array_OBS := (');
for i = 1:nx+nd+np+nu+ny
    fprintf(fout,'%d,',preshiftPred(i)-min(preshiftPred));
end
fprintf(fout,'%d);\n',preshiftPred(nx+nd+np+nu+ny+1)-min(preshiftPred));

fprintf(fout,'\n');

fprintf(fout,'constant N_BIT_MUL_OBS : integer := %d;\n',max(nbit,nbit_coeff));

fprintf(fout,'type mul_in_matrix_OBS is array(0 to nx_OBS+nd_OBS-1) of signed(N_BIT_MUL_OBS-1 downto 0);\n');
fprintf(fout,'type mul_out_matrix_OBS is array(0 to nx_OBS+nd_OBS-1) of signed(2*N_BIT_MUL_OBS-1 downto 0);\n');

fprintf(fout,'\n');

fprintf(fout,'end observer_package;\n');

fclose(fout);

disp(['Generated file ',folder,'observer_package.vhdl']);