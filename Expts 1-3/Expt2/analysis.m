%% analyse results - Online study

clear
clc
close all

load data.mat

% Subjects
subjects = 1:length(data);
%subjects = 1:5;
nratings = 4;
trials_per_task = 100;

% Add helper functions
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})-length(parts{end-1})-1));
addpath('helperFunctions','Mfunctions');
cd(root)
%%
%Flags
compute_measures = 0;
save_measures = 0;

if compute_measures
    for subject = subjects
        clear stim resp correct conf rt1 rt2 task
        
        totalTrials = 0;
        for task_type = 1:2 %1:high contrast, 2:low contrast
            
            task_name = {'high','low'};
            
            trials = totalTrials + [1:trials_per_task];
            totalTrials = totalTrials + trials_per_task;
            
            %Behavior
            stim(trials) = eval(char(strcat('data{1,subject}.stimID_', task_name(task_type), '_contrast')));
            resp(trials) = eval(char(strcat('data{1,subject}.response_', task_name(task_type), '_contrast')));
            correct(trials) = eval(char(strcat('data{1,subject}.correct_', task_name(task_type), '_contrast')));
            conf(trials) = eval(char(strcat('data{1,subject}.confidence_', task_name(task_type), '_contrast')));
            
            %Task parameters
            task(trials) = task_type;
            
            
            %Transform to binary confidence ratings by varying confidence criterion
            %location (Between 1&2 (low), 2&3(med), 3&4(high))
            for criterion = 1:nratings-1
                conf_binary = (conf>criterion)+1;
                
                %Mratio and meta-d'
                meta_MLE = type2_SDT_MLE(stim(task==task_type),resp(task==task_type),conf_binary(task==task_type),nratings, [], 1);
                mratio(subject,criterion, task_type) = meta_MLE.M_ratio;
                meta_d(subject,criterion,task_type) = meta_MLE.meta_da;
                              
                %Phi
                R = corrcoef(correct(task==task_type),conf_binary(task==task_type));
                phi(subject,criterion,task_type) = R(2);
                %diff(subject,criterion) = mean(conf_binary(correct==1,criterion))-mean(conf_binary(correct==0,criterion));
                
                %Type-2 AUC
                [nR_S1, nR_S2] = trials2counts(stim(task==task_type),resp(task==task_type),conf_binary(task==task_type),nratings,[]);
                [t2AUC(subject,criterion,task_type)] = type2ag(nR_S1, nR_S2, 1);
                
                %Type-2 d'
                t2_dprime(subject,criterion,task_type) = data_analysis_resp(correct(task==task_type),conf_binary(task==task_type)-1);
                %t2_dprime(t2_dprime<0)=0;
                
                %d'
                dprime(subject,criterion,task_type) = data_analysis_resp(stim(task==task_type),resp(task==task_type));
                
            end
        end
    end
else
    load metaMeasures_expt2.mat
end


if save_measures
    save('metaMeasures_expt2','phi','t2AUC','mratio','meta_d','t2_dprime');
end


%% Plot

measures = {'phi','t2AUC','meta_d','mratio'};
names = {'\phi','Type-2 AUC','meta-d''','M_{ratio}'};

numMeasures = length(names);

for task_type = 1:2
    figure;
    
    for measure = 1:numMeasures
        y_ = eval(measures{measure});
        y = y_(:,:,task_type);
        
        %Remove outliers
        y_sd = std(reshape(y,[],1));
        y_mean = mean(reshape(y,[],1));
        numOutliers(measure,task_type) = length(find(y>y_mean+3*y_sd | y<y_mean-3*y_sd))
        y(y>y_mean+3*y_sd | y<y_mean-3*y_sd) = NaN;        
        
        subplot(2,2,measure); hold on;
        plot([1 2 3],y,'o')
        
%         axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',6); hold on;
%         subplot(2,2, measure); hold on;
%         [pval] = pairwiseComparisons(y);
%         sem = nanstd(y)./sqrt(length(subjects));
%         barPlotData(y,names{measure},pval,sem,1)
%         
%         if task_type == 1
%             suptitle('Fine discrimination task')
%         else
%             suptitle('Coarse discrimination task')
%         end
%         
    end
end



