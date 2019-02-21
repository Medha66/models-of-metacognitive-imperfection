function [ output_args ] = plot_linearFit(y,ylab,sem)
%Do a linear regression on the metacognition measure with confCriteria as a
%predictor - plot the data with error bars and the fitted line. Indicate
%the value of regression coefficient and R-squared

%columns of y correspond to different conditions

%Other nputs
%ylab = Y axis label, ylimits - limits of the y-axis
%sem = standard error

%Plot data
nCriteria = size(y,2);
x = [1:nCriteria];
h1 = errorbar(x,nanmean(y),sem,'.-','LineWidth',.5,'Markersize',10,'color',[.5 .5 .5]); hold on;
h1 = plot(x,nanmean(y),'.','LineWidth',1.5,'Markersize',15,'color',[0 0 0]);
ylabel({'';ylab;''})
xlabel({'Criterion Location';''},'FontSize',10)
% yl = ylim;
% yl(1) = 0;

%Linear regression
[b,~,~,~,stats] = regress(nanmean(y)',[ones([nCriteria,1]),x'])

%Linear fit results
%h1 = plot(x,x*b(2)+b(1),'-r','linewidth',1.5)
str = strcat('slope = ',num2str(round(b(2),2)),',  p = ',num2str(round(stats(3),2)));

l = legend([h1],str);
set(l,'fontSize',8,'location','best')


end

