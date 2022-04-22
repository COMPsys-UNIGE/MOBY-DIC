function generatePackage(object,folder,FintPred, GintPred, alphaPred, betaPred, nbit, nint, nbit_coeff, nint_coeffPred,initialState,sampling_latency)

ndyn = object.nDyn;

ndec = nbit-nint;

ndec_coeffPred = cell(ndyn,1);
shiftMulPred = cell(ndyn,1);
preshiftPred = cell(ndyn,1);

MM = cell(ndyn,1);

for i=1:ndyn
    ndec_coeffPred{i} = nbit_coeff-nint_coeffPred{i};
    
    
    shiftMulPred{i} = [ndec 0]+ndec_coeffPred{i};
    preshiftPred{i} = max(shiftMulPred{i})-shiftMulPred{i};
    
    MM{i} = [FintPred{i} GintPred{i}];
    
end

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;



fout = fopen([folder,'observer_package.vhdl'],'w');

fprintf(fout,'library IEEE;\n');
fprintf(fout,'use IEEE.STD_LOGIC_1164.all;\n\n');
fprintf(fout,'use IEEE.NUMERIC_STD.all;\n\n');
fprintf(fout,'package observer_package is\n\n');

fprintf(fout,'constant N_BIT_OBS : integer := %d;\n',nbit);
fprintf(fout,'constant N_BIT_COEFF_OBS : integer := %d;\n\n',nbit_coeff);

fprintf(fout,'constant sampling_latency_OBS : integer := %d;\n\n',sampling_latency);

fprintf(fout,'constant nx_OBS : integer := %d;\n',object.getNumberOfStates);
fprintf(fout,'constant np_OBS : integer := %d;\n',object.getNumberOfParameters);
fprintf(fout,'constant nd_OBS : integer := %d;\n',object.getNumberOfUnmeasurableInputs);
fprintf(fout,'constant nu_OBS : integer := %d;\n',object.getNumberOfInputs);
fprintf(fout,'constant ny_OBS : integer := %d;\n\n',object.getNumberOfOutputs);
fprintf(fout,'constant nDyn_OBS : integer := %d;\n\n',object.nDyn);

tmp = cell2mat(preshiftPred);
fprintf(fout,'constant minShiftPred_OBS : integer := %d;\n\n',min(tmp(:)));
fprintf(fout,'constant maxShiftPred_OBS : integer := %d;\n\n',max(tmp(:)));

uno = decimal2signed(1,nbit,0);
nBusMacPredict = nbit+nbit_coeff+ceil(log2(nx+np+nd+nu+ny))+max(tmp(:))-min(tmp(:));
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
fprintf(fout,'type BigObserverMatrixType_OBS is array (0 to nDyn_OBS-1) of observerMatrixType_OBS;\n\n');

alphaPredBin = cell(ndyn,1);


fprintf(fout,'type alphaPredType_OBS is array (0 to nDyn_OBS-1) of std_logic_vector(N_BIT_COEFF_OBS-1 downto 0);\n');
fprintf(fout,'constant alphaPred_OBS : alphaPredType_OBS := (');
for i=1:ndyn-1
    alphaPredBin{i} = decimal2signed(alphaPred{i},nbit_coeff);
    fprintf(fout,'"%s",',alphaPredBin{i}.bin);
end
alphaPredBin{end} = decimal2signed(alphaPred{end},nbit_coeff);
fprintf(fout,'"%s");',alphaPredBin{end}.bin);

fprintf(fout,'type betaPredType_OBS is array (0 to nDyn_OBS-1) of std_logic_vector(N_BIT_COEFF_OBS+nBusMacPredict_OBS-1 downto 0);\n');
fprintf(fout,'constant betaPred_OBS : betaPredType_OBS := (\n\n');

for i=1:ndyn-1
    betaPredBin = decimal2signed(betaPred{i},nBusMacPredict+nbit_coeff,max(shiftMulPred{i})+alphaPredBin{i}.pointposition);
    fprintf(fout,'"%s",',betaPredBin.bin);
end
betaPredBin = decimal2signed(betaPred{end},nBusMacPredict+nbit_coeff,max(shiftMulPred{end})+alphaPredBin{end}.pointposition);
fprintf(fout,'"%s");',betaPredBin.bin);


% Write M matrix
fprintf(fout,'constant predictMatrix : BigObserverMatrixType_OBS := ((( ');
for k=1:ndyn
    M = MM{k};
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
            fprintf(fout,')');
        end
    end
    if k ~= ndyn
        fprintf(fout,',\n((');
    else
        fprintf(fout,');\n\n');
    end
end

fprintf(fout,'\ttype int_shift_array_OBS is array(0 to nx_OBS+nd_OBS+nu_OBS+np_OBS+ny_OBS) of integer;\n');
fprintf(fout,'\ttype dyn_int_shift_matrix_OBS is array(0 to nDyn_OBS-1) of int_shift_array_OBS;\n');

fprintf(fout,'\ttype dyn_int_array_OBS is array(0 to nDyn_OBS-1) of integer;\n');

fprintf(fout,'\tconstant nDecOutPred_OBS : dyn_int_array_OBS := (');
for i=1:ndyn-1
    fprintf(fout,'\t%d,',alphaPredBin{i}.pointposition+max(shiftMulPred{i}));
end
fprintf(fout,'\t%d);\n\n',alphaPredBin{end}.pointposition+max(shiftMulPred{end}));


fprintf(fout,'\n');
fprintf(fout,'\tconstant preShiftPredict_OBS : dyn_int_shift_matrix_OBS := ((');
for k = 1:ndyn
    preshiftPred_k = preshiftPred{k};
    for i = 1:nx+nd+np+nu+ny
        fprintf(fout,'%d,',preshiftPred_k(i)-min(preshiftPred_k));
    end
    fprintf(fout,'%d)\n',preshiftPred_k(nx+nd+np+nu+ny+1)-min(preshiftPred_k));
    if k ~= ndyn
        fprintf(fout,',\n(');
    else
        fprintf(fout,');\n\n');
    end
end
fprintf(fout,'\n');

fprintf(fout,'constant N_BIT_MUL_OBS : integer := %d;\n',max(nbit,nbit_coeff));

fprintf(fout,'type mul_in_matrix_OBS is array(0 to nx_OBS+nd_OBS-1) of signed(N_BIT_MUL_OBS-1 downto 0);\n');
fprintf(fout,'type mul_out_matrix_OBS is array(0 to nx_OBS+nd_OBS-1) of signed(2*N_BIT_MUL_OBS-1 downto 0);\n');

fprintf(fout,'type mul_in_bigmatrix_OBS is array(0 to nx_OBS+np_OBS+nd_OBS-1) of signed(N_BIT_MUL_OBS-1 downto 0);\n');
fprintf(fout,'type mul_out_bigmatrix_OBS is array(0 to nx_OBS+np_OBS+nd_OBS-1) of signed(2*N_BIT_MUL_OBS-1 downto 0);\n');


fprintf(fout,'\n');

fprintf(fout,'end observer_package;\n');

fclose(fout);

disp(['Generated file ',folder,'observer_package.vhdl']);