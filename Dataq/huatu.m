% 画图
figure(1)

plot(nmi,'-ok','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7]);    % 标记点
% hold on
% plot(nmi,'LineWidth',1.5,'Color',[.7 .7 .7]);%实际曲线

xlabel('\fontname{Times New Roman}\fontsize{14} Generations');
ylabel('\fontname{Times New Roman}\fontsize{14} NMI');
% saveas(gcf,'NMI.fig');
% saveas(gcf,'NMI.emf');
% saveas(gcf,'NMI.svg');