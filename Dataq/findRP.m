clc
clear
X = [];
 for i = 1:5
    savefile =  sprintf('CommunityNSGAII_NREG_M2_D1156_%d.mat',i);
    load(savefile)
    for j = 1:length(Population)
        X = [X;Population(j).obj];
    end
 end
for k = 1:4
 for i = 1:5
    savefile =  sprintf('CommunityNSGAII%d_NREG_M2_D1156_%d.mat',k,i);
    load(savefile)
    for j = 1:length(Population)
        X = [X;Population(j).obj];
    end
 end
end
referenceP = [max(X(:,1)),max(X(:,2))]