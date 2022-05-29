% load all dataset and find the maximum of PF
% load real network
clc,clear
str = {'karate','polbooks','football','lesmis','dolphins','karate','polbooks','football','lesmis','dolphins'};
% load('karate.txt'); % EG1 EG6 RN1 RN6
% load('polbooks.txt'); % EG2 EG7 RN2 RN7
% load('football.txt'); % EG3 EG8 RN3 RN8
% load('lesmis.txt'); % EG4 EG9 RN4 RN9
% load('dolphins.txt'); % EG5 EG10 RN5 RN10
%  for i = 1:5
% %     savefile =  sprintf('CommunityNSGAII_NREG2_M2_D1156_%d.mat',i);
%     savefile =  sprintf('NSGAII_NREG1_M2_D1156_%d.mat',i);
%     load(savefile)
%     Score(1,i) = cellfun(@(S)Metric(S,Population,[20,200]),{@HV},'UniformOutput',false);
%     % compute MCC and obtain the best value of MCC
%     MCC(1,i) = measureNR(Population,karate);
%  end
D = [1156,11025,13225,5929,3844,1156,11025,13225,5929,3844];
for k = 10
    for i = 1:3
%         savefile =  sprintf('CommunitySPEA2_NRRN%d_M2_D%d_%d.mat',k,D(k),i);
        savefile =  sprintf('CommunitySPEA2_NREG%d_M2_D%d_%d.mat',k,D(k),i);
%         savefile =  sprintf('CommunityNSGAII_NREG%d_M2_D%d_%d.mat',k,D(k),i);
%             savefile =  sprintf('CommunityNSGAII_NRRN%d_M2_D%d_%d.mat',k,D(k),i);
        load(savefile)
        realnet=load([str{k},'.txt']);
        Score1(k,i) = cellfun(@(S)Metric(S,Population,[20,200]),{@HV},'UniformOutput',false);
        % compute MCC and obtain the best value of MCC
        MCC1(k,i) = measureNR(Population,realnet);
        % compute NMI and obtain the best value of NMI
        NMI1(k,i) = max(nmi);
    end
end
for k = 10
    for i = 1:3
%         savefile =  sprintf('SPEA2_NRRN%d_M2_D%d_%d.mat',k,D(k),i);
        savefile =  sprintf('SPEA2_NREG%d_M2_D%d_%d.mat',k,D(k),i);
%         savefile =  sprintf('NSGAII_NREG%d_M2_D%d_%d.mat',k,D(k),i);
%             savefile =  sprintf('NSGAII_NRRN%d_M2_D%d_%d.mat',k,D(k),i);
        load(savefile)
        realnet=load([str{k},'.txt']);
        Score2(k,i) = cellfun(@(S)Metric(S,Population,[20,200]),{@HV},'UniformOutput',false);
        % compute MCC and obtain the best value of MCC
        MCC2(k,i) = measureNR(Population,realnet);
        % compute NMI and obtain the best value of NMI
        NMI2(k,i) = max(nmi);
    end
end
fprintf('Our proposed MCC:%1.2e(%1.2e)\n',mean(MCC1(end,:)),std(MCC1(end,:)));
fprintf('NSGAII MCC:%1.2e(%1.2e)\n',mean(MCC2(end,:)),std(MCC2(end,:)));
fprintf('Our proposed NMI:%1.2e(%1.2e)\n',mean(NMI1(end,:)),std(NMI1(end,:)));
fprintf('NSGAII NMI:%1.2e(%1.2e)\n',mean(NMI2(end,:)),std(NMI2(end,:)));
disp(['Our proposed MCC:', num2str(MCC1(end,:))]);
disp(['NSGAII MCC:', num2str(MCC2(end,:))]);
disp(['Our proposed NMI:', num2str(NMI1(end,:))]);
disp(['NSGAII NMI:', num2str(NMI2(end,:))]);
% a = cell2mat(Score);
% boxplot(a');
% ylabel('HV','FontSize',14,'FontWeight','bold');
% set(gca,'FontSize',14)
% xlabel('\itq','FontSize',14,'FontWeight','bold');
% set(gca,'XTicklabel',{'1','2','3','5','8'});
% figure;boxplot(MCC');
% ylabel('AUC','FontSize',14,'FontWeight','bold');
% set(gca,'FontSize',14)
% xlabel('\itq','FontSize',14,'FontWeight','bold');
% set(gca,'XTicklabel',{'1','2','3','5','8'});
% set(gca,'XTicklabel',{'2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','2^0','2^1','2^2','2^3','2^4','2^5'});
% % compute HV and obtain median HV and IQRs
% median()
% igr()

function value = Metric(metric,Population,PF)
       % Calculate the metric value of the population
       Feasible     = find(all(Population.cons<=0,2));
       NonDominated = NDSort(Population(Feasible).objs,1) == 1;
       try
          value = metric(Population(Feasible(NonDominated)).objs,PF);
       catch
          value = NaN;
       end
end