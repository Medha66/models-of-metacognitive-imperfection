% Test dependence of mratio on criterion location

clear
clc
close all

% Subjects
subject_id = [1 2 3 5 6 7 9 10 11 15 16 17 19 20 22 23 24 26 27];
subjects = 1:length(subject_id);
nratings = 4;

% Add helper functions
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})-length(parts{end-1})-1));
addpath('helperFunctions','Mfunctions');
cd(root)


%% Add data
addpath('Data')

%Flags
compute_measures = 1;
save_measures = 1;

%Get the counts
if compute_measures
    for subject = 12
        
        %% Load the data
        clear stim resp correct conf rt
        
        % Load the data
        file_name = ['Data/results_s' num2str(subject_id(subject))];
        eval(['load ' file_name '']);
        
        % Loop over all blocks except TMS test blocks
        main_blocks = [4:12,16:24];
        block_length = length(p.data{main_blocks(1)}.stimulus);
        num_blocks = length(main_blocks);
        
        totalTrials = 0;
        for block=1:length(main_blocks)
            trials = totalTrials + [1:block_length];
            totalTrials = totalTrials + block_length;
            stim(trials) = p.data{main_blocks(block)}.stimulus-1; %1: left, 2: right
            resp(trials) = p.data{main_blocks(block)}.response-1;%1: left, 2: right
            correct(trials) = p.data{main_blocks(block)}.correct; %0: error, 1: correct
            conf(trials) = p.data{main_blocks(block)}.confidence; %1-4
        end
        
        
        %Transform to binary confidence ratings by varying confidence criterion
        %location (Between 1&2 (low), 2&3(med), 3&4(high))
        for criterion = 1:nratings-1
            conf_binary(conf>criterion,criterion)=2;
            conf_binary(conf<=criterion,criterion)=1;
            
            %Mratio and meta-d'
            meta_MLE = type2_SDT_MLE(stim,resp,conf_binary(:,criterion)',nratings, [], 1);
            mratio(subject,criterion) = meta_MLE.M_ratio;
            meta_d(subject,criterion) = meta_MLE.meta_da;
            
            %Phi
            R = corrcoef(correct,conf_binary(:,criterion)');
            phi(subject,criterion) = R(2);
            
            
            %Type-2 AUC
            [nR_S1, nR_S2] = trials2counts(stim,resp,conf_binary(:,criterion)',nratings,[]);
            [t2AUC(subject,criterion)] = type2ag(nR_S1, nR_S2, 1);
            
            %Type-2 d'
            t2_dprime(subject,criterion) = data_analysis_resp(correct,conf_binary(:,criterion)'-1);
            t2_dprime(t2_dprime<0)=0;
            
            %d'
            dprime(subject,criterion) = data_analysis_resp(stim,resp);
            
            
        end
        
    end
else
    load('metaMeasures_expt1')
    
end
%%
if save_measures
    save('metaMeasures_expt1','phi','t2AUC','mratio','meta_d','t2_dprime');
end

%% Plot metacognition measures for each criterion location

measures = {'mratio','meta_d','t2AUC','phi'};
names = {'M_{ratio}','meta-d''','Type-2 AUC','\phi',};

numMeasures = length(names);

for measure = 1:numMeasures
    y = eval(measures{measure});
    
    %Remove outliers
    y_sd = std(reshape(y,[],1));
    y_mean = mean(reshape(y,[],1));
    numOutliers(measure) = length(find(y>y_mean+3*y_sd | y<y_mean-3*y_sd))
    y(y>y_mean+3*y_sd | y<y_mean-3*y_sd) = NaN;
    
    subplot(2,2,measure); hold on;
    plot([1 2 3],y,'o')
    
    
%     subplot(2,2, measure); hold on;
%     [pval] = pairwiseComparisons(y);
%     sem = nanstd(y)./sqrt(length(subjects));
%     barPlotData(y,names{measure},pval,sem,1)
end

%%



