function[Hmem,Kmem,Fmem,Gmem,stateMem,FGi_size] = reduceForCircuit(object,H,K,F,G,state,idx)

Hmem = [];
Kmem = [];
Fmem = [];
Gmem = [];
FGi_size = [];

Ftmp = cell(numel(idx),1);
Gtmp = cell(numel(idx),1);
% stateMem = [];
allAddr = [];
edgeWritten = [];

functionWritten = cell(numel(idx),1);

for i=1:numel(state)
    stateMem(i) = state(i);
    
    if ~stateMem(i).leaf
        ii = find(edgeWritten == state(i).addrHK);
        if isempty(ii)
            edgeWritten = [edgeWritten;state(i).addrHK];
            Hmem = [Hmem;H(state(i).addrHK,:)];
            Kmem = [Kmem;K(state(i).addrHK)];
            stateMem(i).addrHK = numel(edgeWritten);
        else
            stateMem(i).addrHK = ii;
        end
    else
        addr = [];
        for j = 1:numel(idx)
            ii = find(functionWritten{j} == state(i).addrFG(idx(j)));
            if isempty(ii)
                functionWritten{j} = [functionWritten{j};state(i).addrFG(idx(j))];
                Ftmp{j} = [Ftmp{j};F(state(i).addrFG(idx(j)),:)];
                Gtmp{j} = [Gtmp{j};G(state(i).addrFG(idx(j)))];
                addr = [addr,numel(functionWritten{j})];
            else
                addr = [addr,ii];
            end
        end
        allAddr = [allAddr;addr];
    end
end    
HKsize = size(Hmem,1);

sumVect = [HKsize];
for i=1:numel(idx)-1
    sumVect = [sumVect,sumVect(i)+size(Ftmp{i},1)];
end


sumVect = repmat(sumVect,size(allAddr,1),1);
allAddr = allAddr + sumVect;

j = 1;
for i=1:numel(stateMem)
     if stateMem(i).leaf
         stateMem(i).addrFG = allAddr(j,:);
         j = j+1;
     end
end
    
for i=1:numel(idx)
    FGi_size = [FGi_size, size(Ftmp{i},1)];
    Fmem = [Fmem ; Ftmp{i}];
    Gmem = [Gmem; Gtmp{i}];
end

