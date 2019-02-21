%% Plot the metaNoise estimates for different criterion values

clear
clc
close all

% Subjects
subjects = 1:20;
nratings = 6;

%Add helper functions
addpath('helperFunctions','fittingResults/')

% Load model fits
load(['fittingResults/fittingResults_normal'])
metaNoise_normal = modelFit.params;
load(['fittingResults/fittingResults_lognormal'])
metaNoise_lognormal = modelFit.params;

%% plot metaNoise against criterion location (pool all contrasts)
measures = {'metaNoise_lognormal','metaNoise_normal'};
names = {'Meta noise parameter (\sigma_{meta})'};
titles = {'Lognormal meta noise','Gaussian meta noise'};


numMeasures = length(names);

for measure = 1:2
    y_ = eval(measures{measure});
    
    
    for criterion = 1:5
        y(:,criterion) = mean(y_(:,:,criterion),2);
    end
    
    %Plot data
    subplot(1,2, measure); hold on;
    [tstat,meanDiff,pval] = pairwiseComparisons(y);
    sem = nanstd(y)./sqrt(length(subjects));
    barPlotData(y,names{1},pval,sem)
    title(titles{measure})
    
    %ANOVA results
    [p,t,stats] = onewayRepmeasuresANOVA(y)
    if p(1) > .01
        p(1) = round(p(1),2)
    end
    
    
end

%% Figure 8 - metaNoise vs criterion location (plot each contrast separately)

colors = {[.8 .1 .1],[.4 .6 .4],[.7 .7 1]};

conf_binSize = .5/6;
for rating = 1:5
    criterion_vect(rating) = .5+(rating)*conf_binSize;    
end

for measure = 1:2
    
    y_ = eval(measures{measure});
    
    %plot each contrast separately
    for condition = 1:3
        for criterion = 1:5
            y(:,criterion) = y_(:,condition,criterion);
        end
        y_sem = nanstd(y)./sqrt(20)
        
        %Plot data
        subplot(2,2, measure); hold on;
        %plot(criterion_vect,mean(y),'o-','linewidth',2,'color',colors{condition}); hold on;
        h(condition) = errorbar(criterion_vect,mean(y), y_sem,'o-','linewidth',2,'color',colors{condition}); hold on;
        %shadedplot(criterion_vect,mean(y)-y_sem,mean(y)+y_sem,colors{condition})
        %alpha(.1)
        ylim([.1,1.1]) 
        set(gca,'XTick',round(criterion_vect,2))
        xlabel('Confidence criterion location','fontsize',16)
        ylabel({'Meta noise parameter'; '(\sigma_{meta})'},'fontsize',16)
        title(titles{measure},'fontsize',18)
        
    end
    
    if measure == 2 
        l = legend(h,'Contrast 1','Contrast 2','Contrast 3')
        set(l,'fontsize',16','orientation','horizontal')
    end
    % 2-way repeated measures ANOVA
    [p,t,stats] = twowayRepmeasuresANOVA(y_)
    
end


