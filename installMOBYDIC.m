disp(' ')
disp('********************************')
disp('* Welcome to MOBY-DIC Toolbox! *')
disp('********************************')
disp(' ')
disp('Checking MPT3 installation...')
disp(' ')
answ = which('mpt_init.m');

abort = false;
if isempty(answ)
    disp('--> ERROR: MPT3 not found!')
    disp('    You can download MPT3 at the following link:')
    disp('    https://www.mpt3.org/')
    disp(' ')
    abort = true;
end

try
    P = Polyhedron(randn(10,2), ones(10,1));
    P.project(5*randn(2,1));
    P.separate(5*randn(2,1));
    P.interiorPoint;
    clear P
catch
    disp('--> ERROR: MPT3 not found or not working properly!')
    disp('    You can download MPT3 at the following link:')
    disp('    https://www.mpt3.org/')
    disp(' ')
    abort = true;
end

if ~abort
    disp('MPT3 ok.')
    disp(' ')
end

disp('Checking Xilinx Vitis HLS installation...')
disp(' ')
evalc('abort = system(''vitis_hls help'')');
if abort ~= 0
    warning('--> ERROR: Xilinx Vitis HLS not found or not working properly!')
    warning('    You cannnot generate VHDL code for implicit MPC controller without Xilinx Vitis HLS')
    warning('    If Vitis HLS is installed, make sure that the installation folder is added to the PATH environment variable')
else
    disp('Xilinx Vitis HLS ok.')
    disp(' ')
end

disp('Adding folders to path...')
disp(' ')

addpath([pwd,'/classes'])
disp([pwd,'/classes'])
addpath([pwd,'/simulink'])
disp([pwd,'/simulink'])
addpath([pwd,'/vhdl'])
disp([pwd,'/vhdl'])
addpath([pwd,'/c'])
disp([pwd,'/c'])
addpath([pwd,'/cpp'])
disp([pwd,'/cpp'])
addpath(genpath([pwd,'/functions']))
disp([pwd,'/functions'])

disp(' ')
disp('Done.')
disp(' ')

status = savepath;

if status == 1
    disp('WARNING: path could not be saved!')
    disp('A startup.m has been generated. Place it in your startup folder.')
    
    fp = fopen('startup.m','w');
    fprintf(fp,'addpath([pwd,''/classes''])\n');
    fprintf(fp,'addpath([pwd,''/simulink''])\n');
    fprintf(fp,'addpath([pwd,''/vhdl''])\n');
    fprintf(fp,'addpath([pwd,''/functions''])\n');
    fclose(fp);
end
disp(' ')

answ = questdlg('Do you want to compile the mex files?', ...
	'Mex files compilation', ...
	'Yes','No','No');

if strcmpi(answ,'yes') 
    
    disp('Compiling mex files...')
    
    cd functions/mex/
    
    mex pwasFunction_alphamex.c -silent
    warning off
    disp('pwasFunction_alphamex.c compiled')
    mex pwasFunction_evalmex.c -silent
    disp('pwasFunction_evalmex.c compiled')
    mex pwagFunction_evalmex.c -silent
    disp('pwagFunction_evalmex.c compiled')
    mex pwagFunction_findRegionmex.c -silent
    disp('pwagFunction_findRegionmex.c compiled')
    mex MPCctrl_evalmex.c -silent
    disp('MPCctrl_evalmex.c compiled')
    mex linear2matrixmex.c -silent
    disp('linear2matrixmex.c compiled')
    warning on
    
    cd ../../
    
    disp(' ')
    disp('Done.')
    disp(' ')
end

disp(' ')
disp('Installation complete.');
answ = questdlg('Do you want to set the toolbox options?', ...
	'Toolbox options', ...
	'Yes','No','No');

if strcmpi(answ,'yes') 
    MOBYDIC_init
end
