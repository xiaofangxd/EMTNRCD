function [MCCB] = measureNR(Population,W)
labels = reshape(W,1,size(W,1)*size(W,1));
for k = 1:length(Population)

[~,~,~,AUC(k)] = perfcurve(labels,Population(k).dec,1);

xp = reshape(Population(k).dec,size(W,1),size(W,1));
T = size(W,1);
for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) <= 0.05
             xp(i,j) = 0;
	    else
			 xp(i,j) = 1;
        end
     end
end
tp = 0;tn = 0;
fn = 0;fp = 0;
% calculate MCC
for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) > 0 && abs(W(i,j)) > 0
           tp = tp+1;
        elseif abs(xp(i,j)) == 0 && abs(W(i,j)) == 0
            tn = tn+1;
        elseif abs(xp(i,j)) > 0 && abs(W(i,j)) == 0
            fp = fp+1;
        else
            fn = fn+1;
        end
    end
end
MCC(k) = (tp*tn-fp*fn)/sqrt((tp+fn)*(tp+fp)*(tn+fp)*(tn+fn));
end
% MCCB = AUC;
MCC1 = max(AUC);
MCCB = MCC1(1);
end