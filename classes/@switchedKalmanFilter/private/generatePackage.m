function generatePackage(object,folder,FintPred, GintPred, alphaPred, betaPred, nbit, nint, nbit_coeff, nint_coeffPred,FintUpdate, GintUpdate, alphaUpdate, betaUpdate, nint_coeffUpdate,initialState,sampling_latency)

ndyn = object.nDyn;

ndec = nbit-nint;

ndec_coeffPred = cell(ndyn,1);
shiftMulPred = cell(ndyn,1);
preshiftPred = cell(ndyn,1);

MMPred = cell(ndyn,1);

ndec_coeffUpdate = cell(ndyn,1);
shiftMulUpdate = cell(ndyn,1);
preshiftUpdate = cell(ndyn,1);

MMUpdate = cell(ndyn,1);

for i=1:ndyn
    ndec_coeffUpdate{i} = nbit_coeff-nint_coeffUpdate{i};
    
    
    shiftMulUpdate{i} = [ndec 0]+ndec_coeffUpdate{i};
    preshiftUpdate{i} = max(shiftMulUpdate{i})-shiftMulUpdate{i};
    
    MMUpdate{i} = [FintUpdate{i} GintUpdate{i}];
    
    ndec_coeffPred{i} = nbit_coeff-nint_coeffPred{i};
    
    
    shiftMulPred{i} = [ndec 0]+ndec_coeffPred{i};
    preshiftPred{i} = max(shiftMulPred{i})-shiftMulPred{i};
    
    MMPred{i} = [FintPred{i} GintPred{i}];
    
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

tmpPredict = cell2mat(preshiftPred);
fprintf(fout,'constant minShiftPred_OBS : integer := %d;\n\n',min(tmpPredict(:)));
fprintf(fout,'constant maxShiftPred_OBS : integer := %d;\n\n',max(tmpPredict(:)));

tmpUpdate = cell2mat(preshiftUpdate);
fprintf(fout,'constant minShiftUpdate_OBS : integer := %d;\n\n',min(tmpUpdate(:)));
fprintf(fout,'constant maxShiftUpdate_OBS : integer := %d;\n\n',max(tmpUpdate(:)));

uno = decimal2signed(1,nbit,0);
nBusMacPredict = nbit+nbit_coeff+ceil(log2(nx+np+nd+nu+ny))+max(tmpPredict(:))-min(tmpPredict(:));
nBusMacUpdate = nbit+nbit_coeff+ceil(log2(nx+np+nd+nu+ny))+max(tmpUpdate(:))-min(tmpUpdate(:));
fprintf(fout,'constant nBusMacPredict_OBS : integer := %d;\n',nBusMacPredict);%N_BIT_OBS+N_BIT_COEFF_OBS+nx_OBS+np_OBS+nd_OBS+nu_OBS+ny_OBS+maxShiftPred_OBS-minShiftPred_OBS;\n\n');
fprintf(fout,'constant nBusMacUpdate_OBS : integer := %d;\n',nBusMacUpdate);%N_BIT_OBS+N_BIT_COEFF_OBS+nx_OBS+np_OBS+nd_OBS+nu_OBS+ny_OBS+maxShiftPred_OBS-minShiftPred_OBS;\n\n');
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
fprintf(fout,'"%s");\n',alphaPredBin{end}.bin);

fprintf(fout,'type betaPredType_OBS is array (0 to nDyn_OBS-1) of std_logic_vector(N_BIT_COEFF_OBS+nBusMacPredict_OBS-1 downto 0);\n');
fprintf(fout,'constant betaPred_OBS : betaPredType_OBS := (');

for i=1:ndyn-1
    betaPredBin = decimal2signed(betaPred{i},nBusMacPredict+nbit_coeff,max(shiftMulPred{i})+alphaPredBin{i}.pointposition);
    fprintf(fout,'"%s",',betaPredBin.bin);
end
betaPredBin = decimal2signed(betaPred{end},nBusMacPredict+nbit_coeff,max(shiftMulPred{end})+alphaPredBin{end}.pointposition);
fprintf(fout,'"%s");\n\n',betaPredBin.bin);

alphaUpdateBin = cell(ndyn,1);

fprintf(fout,'type alphaUpdateType_OBS is array (0 to nDyn_OBS-1) of std_logic_vector(N_BIT_COEFF_OBS-1 downto 0);\n');
fprintf(fout,'constant alphaUpdate_OBS : alphaUpdateType_OBS := (');
for i=1:ndyn-1
    alphaUpdateBin{i} = decimal2signed(alphaUpdate{i},nbit_coeff);
    fprintf(fout,'"%s",',alphaPredBin{i}.bin);
end
alphaUpdateBin{end} = decimal2signed(alphaUpdate{end},nbit_coeff);
fprintf(fout,'"%s");\n',alphaUpdateBin{end}.bin);

fprintf(fout,'type betaUpdateType_OBS is array (0 to nDyn_OBS-1) of std_logic_vector(N_BIT_COEFF_OBS+nBusMacUpdate_OBS-1 downto 0);\n');
fprintf(fout,'constant betaUpdate_OBS : betaUpdateType_OBS := (');

for i=1:ndyn-1
    betaUpdateBin = decimal2signed(betaUpdate{i},nBusMacUpdate+nbit_coeff,max(shiftMulUpdate{i})+alphaUpdateBin{i}.pointposition);
    fprintf(fout,'"%s",',betaUpdateBin.bin);
end
betaUpdateBin = decimal2signed(betaUpdate{end},nBusMacUpdate+nbit_coeff,max(shiftMulUpdate{end})+alphaUpdateBin{end}.pointposition);
fprintf(fout,'"%s");\n\n',betaUpdateBin.bin);


% Write M matrix
fprintf(fout,'constant predictMatrix : BigObserverMatrixType_OBS := ((( ');
for k=1:ndyn
    M = MMPred{k};
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

% Write M matrix
fprintf(fout,'constant updateMatrix : BigObserverMatrixType_OBS := ((( ');
for k=1:ndyn
    M = MMUpdate{k};
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

fprintf(fout,'\tconstant nDecOutUpdate_OBS : dyn_int_array_OBS := (');
for i=1:ndyn-1
    fprintf(fout,'\t%d,',alphaUpdateBin{i}.pointposition+max(shiftMulUpdate{i}));
end
fprintf(fout,'\t%d);\n\n',alphaUpdateBin{end}.pointposition+max(shiftMulUpdate{end}));


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

fprintf(fout,'\n');
fprintf(fout,'\tconstant preShiftUpdate_OBS : dyn_int_shift_matrix_OBS := ((');
for k = 1:ndyn
    preshiftUpdate_k = preshiftUpdate{k};
    for i = 1:nx+nd+np+nu+ny
        fprintf(fout,'%d,',preshiftUpdate_k(i)-min(preshiftUpdate_k));
    end
    fprintf(fout,'%d)\n',preshiftUpdate_k(nx+nd+np+nu+ny+1)-min(preshiftUpdate_k));
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