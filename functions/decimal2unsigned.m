function xint = decimal2unsigned(x,nbit,pointposition)

if ~exist('pointposition','var')
    pointposition = [];
end

if any(x < 0)
    error('Only positive numbers can be converted to unsigned');
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
        xint(i).representation = 'unsigned';
        xint(i).dec = 0;
        xint(i).bin = dec2bin(0,nbit);
    else
        
        xcur = xcur+eps;
        nint = ceil(log2(xcur));
        if isempty(pointposition)
            ndec = nbit-nint;
        else
            ndec = pointposition;
        end
        
        xr = xcur*2^ndec;
        xr = round(xr);
        
        maxval = 2^nbit-1;
        
        xint(i).resolution = nbit;
        xint(i).pointposition = ndec;
        xint(i).representation = 'unsigned';
        if xr > maxval
            xint(i).dec = maxval/2^ndec;
            xint(i).bin = dec2bin(maxval,nbit);
        else
            xint(i).dec = xr/2^ndec;
            xint(i).bin = dec2bin(xr,nbit);
        end
    end

end