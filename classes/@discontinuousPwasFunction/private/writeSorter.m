function writeSorter(object,options)

ndim = object.getDomainDimensions();
folder = options.folder;

if ndim > 4
    fin = fopen([getvhdlpath(),'discontinuousPwasFunction/Sorter8.vhd'],'r');
    ndimmax = 8;
elseif ndim > 2
    fin = fopen([getvhdlpath(),'discontinuousPwasFunction/Sorter4.vhd'],'r');
    ndimmax = 4;
else
    fin = fopen([getvhdlpath(),'discontinuousPwasFunction/Sorter2.vhd'],'r');
    ndimmax = 2;
end

fout = fopen([folder,'Sorter.vhd'],'w');

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
        
        for i = 1:ndim
            fprintf(fout,'\tsignal_x%d <= x(%d);\n',i,i-1);
        end
        for i = ndim+1:ndimmax
            fprintf(fout,'\tsignal_x%d <= (others => ''0'');\n',i);
        end
        
        for i = 1:ndim
            fprintf(fout,'\ty(%d) <= signal_y%d;\n',i-1,i);
        end
        
    else
        
        % ... otherwise write the read line in output file
        fprintf(fout,tline);
        fprintf(fout,'\r\n');
        
    end
    
end

fclose(fin);
fclose(fout);

disp(['Generated file ',folder,'Scale.vhd']);
