function nedge = writeFindDynamic(object,circuit_parameters)

nx = object.nx;
np = object.np;
nu = object.np;
ny = object.ny;
nd = object.nd;


nbit = circuit_parameters.inputResolution;

nbit_coeff = circuit_parameters.coeffResolution;

folder = circuit_parameters.folder;

range = circuit_parameters.range;

inputRange = circuit_parameters.inputRange;

outputRange = circuit_parameters.outputRange;

deltaInADC = inputRange.max-inputRange.min;

meanInADC =  ceil((inputRange.max+inputRange.min)/2);

deltaOutADC = outputRange.max-outputRange.min;

meanOutADC =  ceil((inputRange.max+inputRange.min)/2);

deltaRealIn = [range.xmax(:) - range.xmin(:);...
    range.pmax(:) - range.pmin(:);...
    range.dmax(:) - range.dmin(:)]';

deltaInADC = [outputRange.max(1:nx)-outputRange.min(1:nx);...
    inputRange.max(object.nu+1:object.nu+object.np)-inputRange.min(object.nu+1:object.nu+object.np);...
    outputRange.max(nx+1:end)-outputRange.min(nx+1:end)];

meanRealIn = [(range.xmax(:)+range.xmin(:))/2;...
    (range.pmax(:)+range.pmin(:))/2;...
    (range.dmax(:)+range.dmin(:))/2];


mulScale =  deltaRealIn(:)./deltaInADC(:);
mulScale = mulScale(:)';

alphaZ = diag(mulScale);


ff = object.filters;

range = circuit_parameters.range;
Hd = [eye(object.nx+object.np+object.nd);-eye(object.nx+object.np+object.nd)];
Kd = [range.xmax(:); range.pmax(:); range.dmax(:); - range.xmin(:); -range.pmin(:); -range.dmin(:)];

HH = [];
KK = [];

reg = cell(object.nDyn,1);

for i=1:object.nDyn
    reg{i} = Polyhedron([Hd;ff(i).H],[Kd;ff(i).K]);
    HH = [HH; ff(i).H];
    KK = [KK; ff(i).K];
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

res = zeros(numel(Kred),object.nDyn);

for i=1:object.nDyn
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

fin = fopen([getvhdlpath(),'switchedKalmanPredictor/findDynamic.vhd'],'r');
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
            fprintf(fout,'type HKType_row_OBS is array (0 to nx_OBS+np_OBS+nd_OBS) of signed(N_BIT_COEFF_OBS-1 downto 0);\n');
            if size(Mint,1) == 1 % TO DO mettere a posto il fatto che ci possa essere solo un edge
                fprintf(fout,'type HKType_OBS is array (0 to N_EDGE) of HKType_row_OBS;\n');
            else
                fprintf(fout,'type HKType_OBS is array (0 to N_EDGE-1) of HKType_row_OBS;\n');
            end
            fprintf(fout,'type HK_shift_Type_OBS is array (0 to nx_OBS+np_OBS+nd_OBS) of integer;\n');
            fprintf(fout,'constant HK : HKType_OBS :=(('); 
            % TO DO mettere a posto il fatto che ci possa essere solo un edge
            if size(Mint,1) == 1
                Mint = [Mint; Mint];
            end
            for i=1:size(Mint,1)
                for j=1:size(Mint,2)-1
                    numi = decimal2signed(Mint(i,j),nbit_coeff,0);%nintcoeff(j));
                    fprintf(fout,'"%s",',numi.bin);
                end
                numi = decimal2signed(-Mint(i,end),nbit_coeff,0);%nintcoeff(end));
                fprintf(fout,'"%s")',numi.bin);
                if i ~= size(Mint,1)
                    fprintf(fout,',\n(');
                else
                    fprintf(fout,');\n\n');
                end
            end
            maxDecCoeff = max(nbit_coeff-nintcoeff);
            fprintf(fout,'constant maxSchift_findDyn : integer := %d;\n',max(nbit_coeff-nintcoeff)-min(nbit_coeff-nintcoeff));
            fprintf(fout,'constant shift_findDyn : HK_shift_Type_OBS :=(');
            for i=1:numel(nintcoeff)-1
                fprintf(fout,'%d,',maxDecCoeff-nbit_coeff+nintcoeff(i));
            end
            fprintf(fout,'%d);\n',maxDecCoeff-nbit_coeff+nintcoeff(end));
            
        elseif count == 1
            fprintf(fout,'next_dynamic_process : process(result_reg)\n');
            fprintf(fout,'begin\n');
            fprintf(fout,'case result_reg is\n');
            for i=1:object.nDyn
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