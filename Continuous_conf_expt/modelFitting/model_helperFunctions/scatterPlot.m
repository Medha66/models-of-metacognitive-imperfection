function [ output_args ] = scatterPlot(modelData,exptData,measureName,limits)
%Scatter plots of model and experimental data

model_names = {'Standard','Gamma'};

figure
axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',13); hold on;
for model = 1:size(modelData,2)
    subplot(size(modelData,2),1,model)
    plot(exptData,modelData(:,model),'.','markersize',10,'Color',[rand rand rand],'LineWidth',2); hold on;
    plot(limits,limits,'--k','LineWidth',.2)  
    for subject = 1:length(exptData)
        text(exptData(subject)+.01,modelData(subject,model)+.07,num2str(subject),'fontsize',8)
    end
    ylabel(model_names(model),'FontSize',16)
    xlabel('Experimental Data','FontSize',16)
end
suptitle(measureName)

end

