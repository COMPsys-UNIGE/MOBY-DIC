function xint = decimal2signed(x,nbit,pointposition)

if ~exist('pointposition','var')
    pointposition = [];
end

nx = numel(x);

xint = struct('resolution',cell(nx,1),'pointposition',cell(nx,1),...
    'representation',cell(nx,1),'dec',cell(nx,1),'bin',cell(nx,1));

eps = 1e-8;

for i = 1:nx
    xcur = x(i);
    
    if xcur == 0
        xint(i).resolution = nbit;
        if isempty(pointposition)
            xint(i).pointposition = nbit;
        else
            xint(i).pointposition = pointposition;
        end
        xint(i).representation = 'signed';
        xint(i).dec = 0;
        xint(i).bin = dec2bin(0,nbit);
    end
    
    if xcur > 0
        xcur = xcur+eps;
        nint = ceil(log2(xcur))+1;
        if isempty(pointposition)
            ndec = nbit-nint;
        else
            ndec = pointposition;
        end
        
        xr = xcur*2^ndec;
        xr = round(xr);
        
        maxval = 2^(nbit-1)-1;
        
        xint(i).resolution = nbit;
        xint(i).pointposition = ndec;
        xint(i).representation = 'signed';
        if xr > maxval
            xint(i).dec = maxval/2^ndec;
            xint(i).bin = dec2bin(maxval,nbit);
        else
            xint(i).dec = xr/2^ndec;
            xint(i).bin = dec2bin(xr,nbit);
        end
    end
    if xcur < 0
        
        xcur = -xcur;
        nint = ceil(log2(xcur))+1;
        if isempty(pointposition)
            ndec = nbit-nint;
        else
            ndec = pointposition;
        end
        
        xr = xcur*2^ndec;
        xr = round(xr);
        
        maxval = 2^(nbit-1);
        
        xint(i).resolution = nbit;
        xint(i).pointposition = ndec;
        xint(i).representation = 'signed';
        if xr > maxval
            xint(i).dec = maxval/2^ndec;
            xint(i).bin = dec2bin(maxval,nbit);
        else
            xint(i).dec = -xr/2^ndec;
            xint(i).bin = twoscomplement(-xr,nbit);
        end
        
        
    end
    
end


function nbin = twoscomplement(ndec,nbit)

if any(floor(ndec) ~= ndec)
    error('Only integer numbers can be converted to binary form');
end

if ~exist('nbit','var')
    nbit = 0;
end

ndecpos = ndec(ndec >= 0);
ndecneg = ndec(ndec < 0);

nbitreqpos = max(ceil(log2(ndecpos+1)));
nbitreqneg = max(ceil(log2(abs(ndecneg))))+1;

nbitreq = max(nbitreqpos,nbitreqneg);

if nbit == 0
    nbit = nbitreq;
elseif nbit < nbitreq
    error('The number of bits provided is not sufficient to code the numbers');
end

nbin = repmat('0',numel(ndec),nbit);

for i = 1:numel(ndec)
    % If the number is positive, simply call dec2bin
    if ndec(i) >= 0
        nbin(i,:)  = dec2bin(ndec(i),nbit);

        % If the number is negative use two's complement
    else
        ndec(i) = abs(ndec(i));
        
        % Code the positive number
        nbin(i,:) = dec2bin(ndec(i),nbit);
        
        % Replace 0 with 1 and vicecersa
        for j = 1:size(nbin,2)
            if strcmp(nbin(i,j),'0')
                nbin(i,j) = '1';
            else
                nbin(i,j) = '0';
            end
        end
        
        % Convert to decimal and add one
        ntemp = binary2decimal(nbin(i,:))+1;
        
        % Convert back to binary
        nbin(i,:) = dec2bin(ntemp,nbit);
        
    end
end