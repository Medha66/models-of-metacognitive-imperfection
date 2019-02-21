function [ output_args ] = plotComparisons(modelData,exptData,measureName,ylimits)

%Model data is a matrix - each column corresponds to one model
%Expt data is a vector of the experimental data for that measure

model_names = {'Standard','Gamma (no lapse)','Gamma (lapse)'};

figure
axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',16); hold on;
plot([1:size(modelData,2)], [modelData],'ok')
h1 = bar([mean(modelData)],'LineWidth',.25);
h2 = plot([1:size(modelData,2)], mean(exptData)*ones([1 size(modelData,2)]), '--k','LineWidth',2);
l=legend([h1,h2],'Model Data','Experimental data');
set(l,'FontSize',16)
ylim(ylimits)
alpha(0.2)
ylabel({measureName},'FontSize',20)
xt = get(gca,'XTick');
set(gca,'XTick',1:size(modelData,2)+1,'FontSize',16,'XTickLabel',{model_names{1},model_names{2},model_names{3}})


% for model = 1:size(modelData,2)
%     subplot(2,1,model)
%     h = bar([mean(exptData), mean(modelData(:,model))],'BarWidth',0.5,'LineWidth',1.25); hold on;
%     plot([1 2], [exptData, modelData(:,model)],'ok')
%     ylim(ylimits)
%     alpha(0.2)
%     ylabel(measureName)
%     xt = get(gca,'XTick');
%     set(gca,'XTick',1:3,'FontSize',16,'XTickLabel',{'Expt Data','Model'})
%     title(model_names{model})
% end

end

