function [ output_args ] = barPlotData(y,ylab,pval,sem,printXLabel)

%Used for making a bar plot of the data along with the significance lines
%specified by pval

%columns of y correspond to different conditions

%Other nputs
%ylab = Y axis label, pval - pvalue of the comparisons, ylimits - limits of the y-axis
%sem = standard error

pval(pval>.05)=NaN;
step = 0.2;

for condition = 1:3
    h = bar(condition,nanmean(y(:,condition)),'BarWidth',0.5,'LineWidth',1.25); hold on;
    set(h,'FaceColor',[step step step]);
    step = step+.2;
end
alpha(0.6)
xt = get(gca,'XTick');
set(gca,'XTick',1:3,'FontSize',12)

errorbar(1:3,nanmean(y),sem,'.k','LineWidth',1.75,'Markersize',10)
ylabel({ylab},'FontSize',18)
if printXLabel
    xlabel({'Criterion Location'},'FontSize',18)
end
xlim([.5,3.5])

ymin = min(reshape(y,[],1)); ymax = max(reshape(y,[],1));
if ymin<0 ymin = 0; end
if ymax>2 ymax = 2; end
ylim([ymin,ymax])

sigstar({[1 2],[1 3],[2,3]},[pval(1), pval(2), pval(3)])

end

