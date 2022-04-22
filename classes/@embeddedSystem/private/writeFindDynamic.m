function nedge = writeFindDynamic(object,circuit_parameters)
nbit = circuit_parameters.inputResolution;

nbit_coeff = circuit_parameters.coeffResolution;



folder = circuit_parameters.folder;

range = object.range;


inputRange = circuit_parameters.inputRange;

deltaRealIn = [range.xmax(:) - range.xmin(:);...
    range.pmax(:) - range.pmin(:);...
    range.dmax(:) - range.dmin(:)]';

deltaInADC = [repmat(2^nbit_coeff,object.nx,1);inputRange.max(1:object.np)-inputRange.min(1:object.np);repmat(2^nbit_coeff,object.nd,1);];

meanRealIn = [(range.xmax(:)+range.xmin(:))/2;...
    (range.pmax(:)+range.pmin(:))/2;...
    (range.dmax(:)+range.dmin(:))/2];


mulScale =  deltaRealIn(:)./deltaInADC(:);
mulScale = mulScale(:)';

alphaZ = diag(mulScale);




range = object.range;
Hd = [eye(object.nx+object.np+object.nd);-eye(object.nx+object.np+object.nd)];
Kd = [range.xmax(:); range.pmax(:); range.dmax(:); - range.xmin(:); -range.pmin(:); -range.dmin(:)];

HH = [];
KK = [];

reg = cell(object.dynSys.getNumberOfDynamics,1);
Htot = object.dynSys.getMatrices('H');
Ktot = object.dynSys.getMatrices('K');
for i=1:object.dynSys.getNumberOfDynamics
    H = Htot{i};
    K = Ktot{i};
    reg{i} = Polyhedron([Hd;H],[Kd;K]);
    HH = [HH; H];
    KK = [KK; K];
end



MOBYDICpars = getMOBYDICpars;
tol = MOBYDICpars.roundtol;


HK = [HH KK];

% remove same edge
HK = unique(HK,'rows');

% Create edge set
edge_set = 1:size(HK,1);

% Remove domain edges from edge set
HKd = [Hd Kd];
HK = tol*round(HK/tol);
HKd = tol*round(HKd/tol);
rem = ismember(HK,HKd,'rows');
edge_set(rem) = [];



% remove same edge in opposite description
i = 1;
stop = 0;
while ~stop
    ii = edge_set(i);
    currentHK = HK(ii,:);
    kk = find(ismember(HK,-currentHK,'rows'));
    jj = find(ismember(edge_set,kk));
    edge_set(jj) = [];
    
    if i == numel(edge_set)
        stop = 1;
    end
    i = i+1;
end

HKred = HK(edge_set,:);
Hred = HKred(:,1:end-1);
Kred = HKred(:,end);

% fill res vector
res = zeros(numel(Kred),object.dynSys.getNumberOfDynamics);

for i=1:object.dynSys.getNumberOfDynamics
    centr = reg{i}.chebyCenter;
    if centr.r > 0
        res(:,i) = (Hred*centr.x <= Kred);
    else
        error('Empty dynamic found! Chech H and K matrix');
    end
end


Kred = Kred - Hred*meanRealIn;
Hred = Hred*alphaZ;

D = [inputRange.max(:)'; inputRange.min(:)'];
[Hint, Kint, alpha, beta, nint, nintcoeff] = quantizeData_noBias(Hred,Kred,D,nbit,nbit_coeff);
nedge = numel(Kint);
if beta ~= 0
    error('Error in software! Please report bug to toolbox team!');
end

fin = fopen([getvhdlpath(),'embeddedSystem/findDynamic.vhd'],'r');
fout = fopen([folder,'findDynamic.vhd'],'w');

count = 0;

while 1
    Mint = [Hint, Kint];
    % Read line from input file
    tline = fgetl(fin);
    
    % If there are no characters, break
    if ~ischar(tline)
        break
    end
    
    % If the MATLABGEN section starts, write the VHDL code ...
    if strcmp(tline,'--- BEGIN MATLABGEN ---')
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        if count == 0
            fprintf(fout,'constant N_EDGE : integer := %d;\n',size(Hint,1));
            fprintf(fout,'type HKType_row_ES is array (0 to nx_ES+np_ES+nd_ES) of signed(N_BIT_COEFF_ES-1 downto 0);\n');
            if size(Mint,1) == 1 % TO DO mettere a posto il fatto che ci possa essere solo un edge
                fprintf(fout,'type HKType_ES is array (0 to N_EDGE) of HKType_row_ES;\n');
            else
                fprintf(fout,'type HKType_ES is array (0 to N_EDGE-1) of HKType_row_ES;\n');
            end
            fprintf(fout,'type HK_shift_Type_ES is array (0 to nx_ES+np_ES+nd_ES) of integer;\n');
            fprintf(fout,'constant HK : HKType_ES :=((');
            % TO DO mettere a posto il fatto che ci possa essere solo un edge
            if size(Mint,1) == 1
                Mint = [Mint; Mint];
            end
            for i=1:size(Mint,1)
                for j=1:size(Mint,2)-1
                    numi = decimal2signed(Mint(i,j),nbit_coeff,0);
                    fprintf(fout,'"%s",',numi.bin);
                end
                numi = decimal2signed(-Mint(i,end),nbit_coeff,0);
                fprintf(fout,'"%s")',numi.bin);
                if i ~= size(Mint,1)
                    fprintf(fout,',\n(');
                else
                    fprintf(fout,');\n\n');
                end
            end
            maxDecCoeff = max(nbit_coeff-nintcoeff);
            fprintf(fout,'constant maxSchift_findDyn : integer := %d;\n',max(nbit_coeff-nintcoeff)-min(nbit_coeff-nintcoeff));
            fprintf(fout,'constant shift_findDyn : HK_shift_Type_ES :=(');
            for i=1:numel(nintcoeff)-1
                fprintf(fout,'%d,',maxDecCoeff-nbit_coeff+nintcoeff(i));
            end
            fprintf(fout,'%d);\n',maxDecCoeff-nbit_coeff+nintcoeff(end));
            
        elseif count == 1
            fprintf(fout,'next_dynamic_process : process(result_reg)\n');
            fprintf(fout,'begin\n');
            fprintf(fout,'case result_reg is\n');
            for i=1:object.dynSys.getNumberOfDynamics
                fprintf(fout,'\twhen "%s" =>\n',flipud(num2str(res(:,i))));
                fprintf(fout,'\t\tdynamicOut_next <= %d;\n',i-1);
            end
            fprintf(fout,'\twhen others =>\n');
            fprintf(fout,'\t\tdynamicOut_next <= 0;\n');
            fprintf(fout,'end case;\n');
            fprintf(fout,'end process;\n');
        end
        count = count+1;
        fprintf(fout,'\r\n');
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'kalmanPredictorInterface.vhdl']);

end