%% Compile results from all 5 experiments

clear
clc
close all

%Add helper functions and Mfunctions
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})));
addpath('helperFunctions','Mfunctions');
cd(root)

%Measures of metacognition
% measures = {'phi','t2AUC','meta_d','mratio','t2_dprime'};
% names = 
measures = {'mratio','meta_d','t2AUC','phi'};
names = {'meta-d''/d''','meta-d''','Type-2 AUC','phi'};
numMeasures = length(names);

%% Experiment 1 - TMS study

cd([root,'/Expt1'])

load('metaMeasures_expt1.mat')
subjects = 1:length(mratio);

figure;

for measure = 1:numMeasures
    y = eval(measures{measure});
    
    %Remove outliers
    y_sd = std(reshape(y,[],1));
    y_mean = mean(reshape(y,[],1));
    numOutliers(measure) = length(find(y>y_mean+3*y_sd | y<y_mean-3*y_sd))
    y(y>y_mean+3*y_sd | y<y_mean-3*y_sd) = NaN;
    
    subplot(5,numMeasures, measure); hold on;
    %Pairwise comparisons
    [pval(measure,:),tstat(measure,:)] = pairwiseComparisons(y);
    sem = nanstd(y)./sqrt(length(subjects));
    barPlotData(y,names{measure},pval(measure,:),sem,0)
    %ANOVA
    [p,t,stats] = onewayRepmeasuresANOVA(y,'Confidence criterion');
    
    %Save F and p values
    F(measure,1) = t{2,6}; P(measure,1) = p(1);
    
end

expt_count = 1;

%% Experiments 2a and 2b - Online study - fine and coarse discrimination task

cd([root,'/Expt2'])

load('metaMeasures_expt2.mat')
subjects = 1:length(mratio);
clear numOutliers

for task_type = 1:2
    clear del_index
    for measure = 1:numMeasures
        y_ = eval(measures{measure});
        y = y_(:,:,task_type);
        
        %Remove outliers
        y_sd = std(reshape(y,[],1));
        y_mean = mean(reshape(y,[],1));
        numOutliers(measure,task_type) = length(find(y>y_mean+3*y_sd | y<y_mean-3*y_sd));
        y(y>y_mean+3*y_sd | y<y_mean-3*y_sd) = NaN;        
        
        count = 1; 
        
        %Delete rows with nans
        for i = 1:length(y)
            a = y(i,:);
            if any(isnan(a))
                del_index(count) = i;
                count = count + 1;
            end
        end
        y(del_index,:) = [];
              
        axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',6); hold on;
        subplot(5,numMeasures, numMeasures*expt_count+measure); hold on;
        %Pairwise comparisons
        [pval(measure,:),tstat(measure,:)] = pairwiseComparisons(y);
        sem = nanstd(y)./sqrt(length(subjects));
        barPlotData(y,names{measure},pval(measure,:),sem,0)
        %ANOVA
        [p,t,stats] = onewayRepmeasuresANOVA(y,'Confidence criterion');
        
        %Save p and F values
        F(measure,1) = t{2,6}; P(measure,1) = p(1);
        
        
    end
    expt_count = expt_count+1;
end



%% Experiments 3a and 3b - Confidence leak study - color and letter identification tasks

cd([root,'/Expt3'])

load('metaMeasures_expt3.mat')
subjects = 1:length(mratio);

for task_type = 1:2
    clear del_index
    for measure = 1:numMeasures
        y_ = eval(measures{measure});
        y = y_(:,:,task_type);
        
        %Remove outliers
        y_sd = std(reshape(y,[],1));
        y_mean = mean(reshape(y,[],1));
        numOutliers(measure,task_type) = length(find(y>y_mean+3*y_sd | y<y_mean-3*y_sd))
        y(y>y_mean+3*y_sd | y<y_mean-3*y_sd) = NaN; 
    
         count = 1; 
        %Delete rows with nans
        for i = 1:length(y)
            a = y(i,:);
            if any(isnan(a))
                del_index(count) = i;
                count = count + 1;
            end
        end
        if count > 1
        y(del_index,:) = []
        end
        
        axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',6); hold on;
        subplot(5,numMeasures, numMeasures*expt_count+measure); hold on;
        %Pairwise comparisons
        [pval(measure,:),tstat(measure,:)] = pairwiseComparisons(y);
        sem = nanstd(y)./sqrt(length(subjects));
        barPlotData(y,names{measure},pval(measure,:),sem,task_type-1)
        
        % ANOVA
        [p,t,stats] = onewayRepmeasuresANOVA(y,'Confidence criterion');
        
        %Save p and F values
        F(measure,1) = t{2,6}; P(measure,1) = p(1);
        
        
    end
    expt_count = expt_count+1;
end





