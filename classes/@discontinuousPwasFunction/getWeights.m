function w = getWeights(object)
%GETWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

w = cell(object.nDyn,1);

for i=1:object.nDyn
    w{i} = object.functions(i).pwasFunction.getWeights;
end


end

