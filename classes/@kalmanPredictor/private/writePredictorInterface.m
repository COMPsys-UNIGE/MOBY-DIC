function writePredictorInterface(object,folder)
nx = object.getNumberOfStates;
nd = object.getNumberOfUnmeasurableInputs;
np = object.getNumberOfParameters;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;


fin = fopen([getvhdlpath(),'kalmanPredictor/kalmanPredictorInterface.vhdl'],'r');
fout = fopen([folder,'kalmanPredictorInterface.vhdl'],'w');

count = 0;

while 1
    
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
            first = 1;
            for i=1:nu
                if ~first
                    fprintf(fout,';\n');
                end
                fprintf(fout,'\tu%d : in std_logic_vector(N_BIT_OBS-1 downto 0)',i);
                first = 0;
            end
            for i=1:np
                if ~first
                    fprintf(fout,';\n');
                end
                fprintf(fout,'\tp%d : in std_logic_vector(N_BIT_OBS-1 downto 0)',i);
                first = 0;
            end
            for i=1:ny
                if ~first
                    fprintf(fout,';\n');
                end
                fprintf(fout,'\ty%d : in std_logic_vector(N_BIT_OBS-1 downto 0)',i);
                first = 0;
            end
            for i=1:nx
                if ~first
                    fprintf(fout,';\n');
                end
                fprintf(fout,'\tx%d_stim : out std_logic_vector(N_BIT_OBS-1 downto 0)',i);
                first = 0;
            end
            for i=1:nd
                if ~first
                    fprintf(fout,';\n');
                end
                fprintf(fout,'\td%d_stim : out std_logic_vector(N_BIT_OBS-1 downto 0)',i);
                first = 0;
            end
        elseif count == 1
            %Assign input
            for i=1:nu
                fprintf(fout,'u_in(%d) <= u%d;\n',i-1,i);
            end
            for i=1:np
                fprintf(fout,'u_in(%d) <= p%d;\n',i-1+nu,i);
            end
            for i=1:ny
                fprintf(fout,'u_in(%d) <= y%d;\n',i-1+nu+np,i);
            end
            fprintf(fout,'\n');
            %Assign output
            for i=1:nx
                fprintf(fout,'x%d_stim <= x_out_reg(%d);\n',i,i-1);
            end
            for i=1:nd
                fprintf(fout,'d%d_stim <= x_out_reg(%d);\n',i,i-1+nx);
            end
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

