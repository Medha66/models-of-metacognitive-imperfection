function [ output_args ] = barPlotData(y,ylab,pval,sem)

%Used for making a bar plot of the data along with the significance lines
%specified by pval

%columns of y correspond to different conditions

%Other nputs
%ylab = Y axis label, pval - pvalue of the comparisons, ylimits - limits of the y-axis
%sem = standard error

pval(pval>.05)=NaN;
step = 0.2;

for condition = 1:5
    h = bar(condition,nanmean(y(:,condition)),'BarWidth',0.5,'LineWidth',1.25); hold on;
    set(h,'FaceColor',[step step step]);
    step = step+.2;
end
alpha(0.8)
xt = get(gca,'XTick');
set(gca,'XTick',1:5,'FontSize',18)

errorbar(1:5,nanmean(y),sem,'.k','LineWidth',1.75,'Markersize',10)
ylabel({ylab},'FontSize',18)
xlabel({'Criterion Location';''},'FontSize',18)
xlim([0,5.5])
% ymax = max(reshape(y,[],1)); ymin = min(reshape(y,[],1));
% ylim([0,ymax+.05*(ymax)])
sigstar({[1 2],[1 3],[1,4],[1,5]},[pval(1), pval(2), pval(3),pval(4),])

end

