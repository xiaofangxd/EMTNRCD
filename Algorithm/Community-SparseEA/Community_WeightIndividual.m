function obj = Community_WeightIndividual(GlobalDummy)
            
if nargin > 0
    xPrimeSize = GlobalDummy.D;
    tempD = sqrt(GlobalDummy.Global.D);
    x = reshape(GlobalDummy.xPrime.dec, tempD, tempD);
    if isempty(GlobalDummy.Index)
        xPrimeVar = GlobalDummy.xPrime.dec;
        xIndex = GlobalDummy.xIndex;
        yIndex = GlobalDummy.yIndex;
        xg = xPrimeVar(xIndex,yIndex);
    else
        Index = GlobalDummy.Index;
        xg = x(Index,Index);
        xg = reshape(xg,1,GlobalDummy.D);
    end
    obj = xg;

    for i = 1:GlobalDummy.N-1
        xgroups = xg;              
        for j = 1:xPrimeSize            
            if rand > 0.8               
                if xgroups(j) == 1                   
                    xgroups(j) = 0;%% can be changed if encoding 'binary'
                else
                    xgroups(j) = 1;
                end
            end
        end
        obj = [obj;xgroups];
    end
end
end