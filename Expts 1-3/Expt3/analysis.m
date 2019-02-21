%% Analyse metacognition measures at different criterion locations

clear
close all

% Select subjects (S9 is completely at chance)
selected_subjects = [1:8,10:27];
%selected_subjects = [1];
number_subjects = length(selected_subjects);
nratings = 4;


% Add path to helper functions
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})-length(parts{end-1})-1));
addpath('helperFunctions','Mfunctions');
cd(root)
addpath('Data');

% Flags
compute_measures = 0;
save_measures = 0;

% Loop through all subjects

if compute_measures
    for subject=1:number_subjects
        
        %Load the data
        file_name = ['Data/results_s' num2str(selected_subjects(subject)) ''];
        eval(['load ' file_name '']);
        
        %Size of blocks and runs
        number_blocks = size(p.rt,1);
        number_trials_per_block = size(p.rt,2);
        trials_per_task = number_blocks*number_trials_per_block;
        
        totalTrials = 0;
        
        for task_type = 1:2 %1 - letter identity, 2 - color
            
            task_name = {'letter','color'};
            
            trials = totalTrials + [1:trials_per_task];
            totalTrials = totalTrials + trials_per_task;
            
            %% Determine the vectors for stimulus, response, and confidence for each task
            if task_type == 1
                stim(trials) = 1 - [p.x_es{1}, p.x_es{2}, p.x_es{3}, p.x_es{4}]'; %0:X, 1:O
                resp(trials) = reshape(p.answers(:,:,1)',400,1) - 1; %0:X, 1:O
                correct(trials) = stim(trials) == resp(trials);
                conf(trials) = reshape(p.answers(:,:,2)',400,1); %letter identity task
                task(trials) = task_type;
                
            else
                stim(trials) = 1 - [p.red{1}, p.red{2}, p.red{3}, p.red{4}]';  %0:red, 1:blue
                resp(trials) = reshape(p.answers(:,:,4)',400,1) - 1; %0:red, 1:blue
                correct(trials) = stim(trials) == resp(trials);
                conf(trials) = reshape(p.answers(:,:,5)',400,1); %color task
                task(trials) = task_type;
            end
            
            
            %Metacognition measures for each criterion location
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
                
                %d'
                dprime(subject,criterion,task_type) = data_analysis_resp(stim(task==task_type),resp(task==task_type));
                
                %Type-2 d'
                t2_dprime(subject,criterion,task_type) = data_analysis_resp(correct(task==task_type),conf_binary(task==task_type)-1);
                %t2_dprime(t2_dprime<0)=0;
                
            end
        end
        
        
    end
else
    load('metaMeasures_expt3.mat')
end
%Remove mratio < 0
%mratio(mratio<0) = NaN;


if save_measures
    save('metaMeasures_expt3','phi','t2AUC','mratio','meta_d','t2_dprime');
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
        
%         subplot(2,2,measure); hold on;
%         plot([1 2 3],y,'o')
        
        axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',6); hold on;
        subplot(2,2, measure); hold on;
        [pval] = pairwiseComparisons(y);
        sem = nanstd(y)./sqrt(number_subjects);
        barPlotData(y,names{measure},pval,sem,1)
        
        if task_type == 1
            suptitle('Letter identitiy task')
        else
            suptitle('Color identity task')
        end
        
    end
end
