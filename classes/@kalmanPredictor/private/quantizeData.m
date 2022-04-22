%% QUANTIZE DATA
% Quantize the coefficients of the PWAG function
%
% DESCRIPTION
%
% This function converts all coefficients of the PWAG function in fixed
% point representation and provides in output all the needed information.
% The quantization procedure is explained below with a simple example.
%
% The operations which must be performed in the circuit are either an edge
% evaluation:
%
% h1 x1 + h2 x2 + k         (1)
%
% or a function computation:
%
% y = f1 x1 + f2 x2 + g     (2)
%
% All terms h1, h2, k, f1, f2, g must be converted in fixed point
% representation with nbit_coeff bits, while for x1 and x2 nbit bits are
% used. All coefficients h1, h2, f1, f2, k, g are stored in matrices H, F,
% K and G.
%
% Data are represented in two's complement coding, therefore, with
% nbit_coeff bits it is possible to represent values between
%
% -2^(nbit_coeff-1) and 2^(nbit_coeff-1)-1      (3).
%
% The first thing to do is to center the values of G around the origin, by
% introducing an offset GG. Expression (2) becomes therefore 
%
% y = f1 x1 + f2 x2 + g - GG    (4)
%
% Secondly F and G are multiplied by a scaling factor alpha such that they 
% completely cover range (3) (i.e. either the maximum or the minimum value
% in F and G coicides with a boundary of the range.
% (4) is then transformed in
%
% y = alpha (f1 x1 + f2 x2 + g - GG)    (5)
%
% When evaluating (1), only the sign is of interest, so it is possible to
% scale each inequality in a different way in order to cover range (3),
% without changing the result of the operation.
%
% From now on, H, K, F and G denote the scaled matrices.
%
% The minimum and maximum values of f1, h1 are now retrieved, in order to
% decide how many bits of integer part are necessary to represent all them.
% The same is done for h2, f2 and k, g. Each component of the computation
% is therefore represented with a possibly different number of decimal
% bits (dhf1, dhf2 and dkg). 
%
% In order to have integer numbers, h1 and f1 are all multiplied by 2^dhf1,
% h2 and f2 by 2^dhf2. Moreover k and g must be multiplied by 2^dkg.
% In order to make computation coherent x1 must be multiplied by
% 2^dkg/2^dhf1 and x2 must be multiplied by 2^dkg/2^dhf2. 
%
% SYNTAX
%
% [Hint Kint Fint Gint A alpha beta nint nintcoeff ndecout] = ...
%   quantizeData(H,K,F,G,D,nbit,nbitcoeff)
%
% H,K,F,G are the coefficient matrices
% D is the domain
% nbit the number of bits to be used to represent inputs
% nbitcoeff the number of bits to be used to represent coefficients
% Hint, Kint, Fint, Gint are the quantized matrices
% A is the scaling factor for the input xcir = A x
% alpha and beta are the scaling factor for the output f = alpha fcir + beta
% nint and nintcoeff are the bits coding the integer part of inputs and
% coefficients
% ndecout is the number of bits to code the decimal part of the output
%
% ACKNOWLEDGEMENTS
%
% % Copyright is with the following author(s):
%
% (C) 2012 Alberto Oliveri, University of Genoa, Italy, alberto.oliveri@unige.it   

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------



function [Fint, Gint, alpha, beta, nint, nintcoeff] = ...
    quantizeData(F,G,D,nbit,nbitcoeff)


% Minimum and maximum value of G
minG = min(G);
maxG = max(G);

% Offset used to center values of G around the origin
offset = (maxG+minG)/2;

% Center values of G around the origin (to have maximum precision in two's
% complement)
G = G-offset;

% Put together matrices F and G
FG = [F G];

% Range representable with nbitcoeff bits in two's complement
rangemin = -2^(nbitcoeff-1);
rangemax = 2^(nbitcoeff-1)-1;

% Take maximum and minimum value of FG
minFG = min(FG(:));
maxFG = max(FG(:));

% Normalize elements of FG in order to completely fill this range
alphaFGmax = abs(rangemax/maxFG);
alphaFGmin = abs(rangemin/minFG);
alphaFG = min(alphaFGmin,alphaFGmax);
% find first power of two less then alphaFG
esp = (ceil(log2(alphaFG+1))-1);
while 1
    if alphaFG >= 2^esp
        break;       
    end
    esp = esp-1;
end
alphaFG = 2^esp; % take the first power of 2 less then alphaG
FG = FG*alphaFG;

% Separate F and G
F = FG(:,1:end-1);
G = FG(:,end);


% Find minimum and maximum values for each column of HF
minF = min(F);
maxF = max(F);

% Find number of bits necessary to represent minimum and maximum values in
% two's complement
nintcoeffmin = ceil(log2(-minF)+1);
nintcoeffmax = ceil(log2(maxF+1)+1);

% Remove complex values (they arise if minHF is positive or maxHF is
% negative)
nintcoeffmin(nintcoeffmin ~= conj(nintcoeffmin)) = NaN;
nintcoeffmax(nintcoeffmax ~= conj(nintcoeffmax)) = NaN;

% Take maximum value
nintcoeffF = max(nintcoeffmin,nintcoeffmax);

% Number of bits remaining for the decimal part
ndeccoeffF = nbitcoeff-nintcoeffF;

% Value used to transform values of F into integers
scaleF = 2.^ndeccoeffF;

% Scale HF
F = F.*repmat(scaleF,size(F,1),1);

% Find minimum and maximum values for each column of KG
minG = min(G);
maxG = max(G);

% Find number of bits necessary to represent minimum and maximum values in
% two's complement
nintcoeffmin = ceil(log2(-minG)+1);
nintcoeffmax = ceil(log2(maxG+1)+1);

% Remove complex values (they arise if minKG is positive or maxKG is
% negative)
nintcoeffmin(nintcoeffmin ~= conj(nintcoeffmin)) = NaN;
nintcoeffmax(nintcoeffmax ~= conj(nintcoeffmax)) = NaN;

% Take maximum value
nintcoeffG = max(nintcoeffmin,nintcoeffmax);

% Number of bits remaining for the decimal part
ndeccoeffG = nbitcoeff-nintcoeffG;

% Value used to transform values of KG into integers
scaleG = 2.^ndeccoeffG;

% Scale KG
G = G*scaleG;


% Round all values
Fint = round(F);
Gint = round(G);

% Saturate values out of range (if any)
Fint(Fint<rangemin) = rangemin;
Fint(Fint>rangemax) = rangemax;
Gint(Gint<rangemin) = rangemin;
Gint(Gint>rangemax) = rangemax;

% % Compute matrices necessary to scale inputs
% % x_cir = A x
% A = scaleG./scaleF;
% 
% % Scale input domain
% DD = D.*repmat(A,2,1);

% Minimum and maximum velues of D
minD = D(1,:);
maxD = D(2,:);

% % Find number of bits necessary to represent minimum and maximum values in
% % two's complement
% nintmin = ceil(log2(-minD)+1);
% nintmax = ceil(log2(maxD+1)+1);
% 
% % Remove complex values (they arise if minD is positive or maxD is
% % negative)
% nintmin(nintmin ~= conj(nintmin)) = NaN;
% nintmax(nintmax ~= conj(nintmax)) = NaN;
% 
% % Take maximum value
% nint = max(nintmin,nintmax);

nint = repmat(nbit,1,size(minD,2));

% if any(nint>nbit)
%     error(['At least ',num2str(max(nint)),' bits are necessary to represent inputs']);
% end

% Number of bits remaining for the decimal part
ndec = nbit-nint;


% Compute matrices necessary to scale output
% f = alpha f_cir + beta
alpha = 1/(alphaFG);
beta = offset;

% A = diag(A);

% Number of bits
nintcoeff = [nintcoeffF nintcoeffG];


% Relative error on F
relerrF = sum(abs(F(:)-Fint(:)))./sum(abs(F(:)))*100;
% Relative error on G
relerrG = sum(abs(G(:)-Gint(:)))./sum(abs(G(:)))*100;

% Maximum quantization error for input
qerr_D = 2.^-ndec;

% Show report
disp(' ')
disp('Relative error on matrix F:')
disp([num2str(relerrF),'%'])
disp('Relative error on matrix G:')
disp([num2str(relerrG),'%'])
disp('Relative quantization error on input points (for each dimension):')
disp(['[',num2str(qerr_D./abs(D(2,:)-D(1,:))*100),'] %'])
disp(' ')


% %%
% 
% % Put together matrices H and K
% HK = [H K];
% 
% % Loop on all rows of HK
% for i = 1:sizeH
%     
%     currrow = HK(i,:);
%     % Take maximum and minimum value of currrow
%     mincurrrow = min(currrow(:));
%     maxcurrrow = max(currrow(:));
%     
%     % Normalize elements of currrow in order to completely fill the range
%     alphacurrrowmax = rangemax/maxcurrrow;
%     alphacurrrowmin = rangemin/mincurrrow;
%     alphacurrrowmin(alphacurrrowmin<0) = Inf;
%     alphacurrrowmax(alphacurrrowmax<0) = Inf;
%     alphacurrrow = min(alphacurrrowmin,alphacurrrowmax);
%     currrow = currrow*alphacurrrow;
%     
%     HK(i,:) = currrow;
%     
% end
% 
% % Separate H and K
% H = HK(:,1:end-1);
% K = HK(:,end);
% 
% % Put together matrices K and G
% KG = [K; G];
% 
% % Find minimum and maximum values for each column of KG
% minKG = min(KG);
% maxKG = max(KG);
% 
% % Find number of bits necessary to represent minimum and maximum values in
% % two's complement
% nintcoeffmin = ceil(log2(-minKG)+1);
% nintcoeffmax = ceil(log2(maxKG)+1);
% 
% % Remove complex values (they arise if minKG is positive or maxKG is
% % negative)
% nintcoeffmin(nintcoeffmin ~= conj(nintcoeffmin)) = NaN;
% nintcoeffmax(nintcoeffmax ~= conj(nintcoeffmax)) = NaN;
% 
% % Take maximum value
% nintcoeffKG = max(nintcoeffmin,nintcoeffmax);
% 
% % Number of bits remaining for the decimal part
% ndeccoeffKG = nbitcoeff-nintcoeffKG;
% 
% % Value used to transform values of KG into integers
% scaleKG = 2.^ndeccoeffKG;
% 
% % Scale KG
% KG = KG*scaleKG;
% 
% % Separate K and G
% K = KG(1:sizeH,:);
% G = KG(sizeH+1:end,:);
% 
% % Put together matrices H and F
% HF = [H; F];
% 
% % Find minimum and maximum values for each column of HF
% minHF = min(HF);
% maxHF = max(HF);
% 
% % Find number of bits necessary to represent minimum and maximum values in
% % two's complement
% nintcoeffmin = ceil(log2(-minHF)+1);
% nintcoeffmax = ceil(log2(maxHF)+1);
% 
% % Remove complex values (they arise if minHF is positive or maxHF is
% % negative)
% nintcoeffmin(nintcoeffmin ~= conj(nintcoeffmin)) = NaN;
% nintcoeffmax(nintcoeffmax ~= conj(nintcoeffmax)) = NaN;
% 
% % Take maximum value
% nintcoeffHF = max(nintcoeffmin,nintcoeffmax);
% 
% % Number of bits remaining for the decimal part
% ndeccoeffHF = nbitcoeff-nintcoeffHF;
% 
% % Value used to transform values of HF into integers
% scaleHF = 2.^ndeccoeffHF;
% 
% % Scale HF
% HF = HF.*repmat(scaleHF,size(HF,1),1);
% 
% % Separate H and F
% H = HF(1:sizeH,:);
% F = HF(sizeH+1:end,:);
% 
% % Round all values
% Hint = round(H);
% Kint = round(K);
% Fint = round(F);
% Gint = round(G);
% 
% % Saturate values out of range (if any)
% Hint(Hint<rangemin) = rangemin;
% Hint(Hint>rangemax) = rangemax;
% Kint(Kint<rangemin) = rangemin;
% Kint(Kint>rangemax) = rangemax;
% Fint(Fint<rangemin) = rangemin;
% Fint(Fint>rangemax) = rangemax;
% Gint(Gint<rangemin) = rangemin;
% Gint(Gint>rangemax) = rangemax;
% 
% % Compute matrices necessary to scale inputs
% % x_cir = A x
% A = scaleKG./scaleHF;
% 
% % Scale input domain
% D = D.*repmat(A,2,1);
% 
% % Minimum and maximum velues of D
% minD = D(1,:);
% maxD = D(2,:);
% 
% % Find number of bits necessary to represent minimum and maximum values in
% % two's complement
% nintmin = ceil(log2(-minD)+1);
% nintmax = ceil(log2(maxD)+1);
% 
% % Remove complex values (they arise if minD is positive or maxD is
% % negative)
% nintmin(nintmin ~= conj(nintmin)) = NaN;
% nintmax(nintmax ~= conj(nintmax)) = NaN;
% 
% % Take maximum value
% nint = max(nintmin,nintmax);
% 
% if any(nint>nbit)
%     error(['At least ',num2str(max(nint)),' bits are necessary to represent inputs']);
% end
% 
% % Number of bits remaining for the decimal part
% ndec = nbit-nint;
% 
% % Number of bits of decimal part for the output
% ndecout = max(ndec);
% 
% % Compute matrices necessary to scale output
% % f = alpha f_cir + beta
% alpha = 1/(alphaFG*scaleKG);
% beta = offset;
% 
% A = diag(A);
% 
% % Number of bits
% nintcoeff = [nintcoeffHF nintcoeffKG];
% 
% % Relative error on H
% relerrH = sum(abs(H(:)-Hint(:)))./sum(abs(H(:)))*100;
% % Relative error on K
% relerrK = sum(abs(K(:)-Kint(:)))./sum(abs(K(:)))*100;
% % Relative error on F
% relerrF = sum(abs(F(:)-Fint(:)))./sum(abs(F(:)))*100;
% % Relative error on G
% relerrG = sum(abs(G(:)-Gint(:)))./sum(abs(G(:)))*100;
% 
% % Maximum quantization error for input
% qerr_D = 2.^-ndec;
% 
% % Show report
% disp(' ')
% disp('Relative error on matrix H:')
% disp([num2str(relerrH),'%'])
% disp('Relative error on matrix K:')
% disp([num2str(relerrK),'%'])
% disp('Relative error on matrix F:')
% disp([num2str(relerrF),'%'])
% disp('Relative error on matrix G:')
% disp([num2str(relerrG),'%'])
% disp('Relative quantization error on input points (for each dimension):')
% disp(['[',num2str(qerr_D./abs(D(2,:)-D(1,:))*100),'] %'])
% disp(' ')