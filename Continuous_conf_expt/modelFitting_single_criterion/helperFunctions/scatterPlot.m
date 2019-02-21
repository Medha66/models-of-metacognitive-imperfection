function [ output_args ] = scatterPlot(exptData, modelData, title, conditionName, modelName)
%Plots the fits from model against exptData
%exptData and modelData - each column represents one condition and each condition is
%plotted as a subplot

numConditions = size(modelData,2);

for condition = 1:numConditions
    
    axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',12); hold on;
    limits = [min(exptData(:,condition)),max(exptData(:,condition))];
    subplot(numConditions,1,condition)
    plot(exptData(:,condition),modelData(:,condition),'.','markersize',15,'Color',[rand rand rand],'LineWidth',2); hold on;
    plot(limits,limits,'--k','LineWidth',.2)
    
    %     for subject = 1:size(exptData,1)
    %         text(exptData(subject,condition)+.01,modelData(subject,condition)+.07,num2str(subject),'fontsize',8)
    %     end
    
    ylabel(['Model Fit'],'FontSize',12)
    xlabel('Experimental Data','FontSize',12)
    l = legend([conditionName ,num2str(condition)]);
    set(l,'location','southeast')
end
suptitle([title,':  ',modelName])


end

