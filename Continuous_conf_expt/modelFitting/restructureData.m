% Restructure data for model fitting

clear
clc
close all

% Subjects
subject_id = [1:13,15,16,18:22];
subjects = 1:length(subject_id);
nratings = 6;
conf_binSize = .5/nratings;
numSessions = 3;

%Flags
compute_measures = 1;
save_measures = 1;

%Add model helper functions
addpath(genpath(fullfile(pwd, 'model_helperFunctions')));

% Add data folder to path
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})));
addpath('Data');
cd(root)

% Add Mfunctions to path
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})-length(parts{end-1})-1));
addpath('Mfunctions');
cd(root)

% Get the counts
for subject = subjects
    clear stim resp correct conf rt contrast sessionID conf_ conf
    
    totalTrials = 0;
    for session = 1:numSessions
        load(['results_sub',num2str(subject_id(subject)),'_s',num2str(session),'.mat'])
        [m,n] = size(p.data);
        numBlocks = m*n;
        main_blocks = [1:numBlocks];
        block_length = length(p.data{1,1}.stim);
        
        for block=1:length(main_blocks)
            trials = totalTrials + [1:block_length];
            totalTrials = totalTrials + block_length;
            stim(trials) = p.data{main_blocks(block)}.stim-1; %1: left, 2: right
            resp(trials) = p.data{main_blocks(block)}.answer-1;%1: left, 2: right
            correct(trials) = p.data{main_blocks(block)}.correct; %0: error, 1: correct
            conf_cont(trials) = p.data{main_blocks(block)}.conf; %.5 - 1
            contrast(trials) = p.data{main_blocks(block)}.contrast; %1-3
            sessionID(trials) = session;
        end
    end    
    
    %Transform conf from 50 point scale to 6 point scale
    for rating = 1:nratings
        conf(conf_cont > .4999+(rating-1)*conf_binSize & conf_cont<= .5+(rating)*conf_binSize) = rating;
    end
    
    %Get counts for each response type
    for condition = 1:3
        for stimulus = 1:2
            for response = 1:2
                for confidence = 1:nratings
                    if response == 1
                        dataCounts(stimulus,7-confidence,subject,condition) = length(find(stim(contrast==condition) == stimulus-1  & resp(contrast==condition) == response-1 & conf(contrast==condition)==confidence));
                    elseif response == 2
                        dataCounts(stimulus,confidence+6,subject,condition) = length(find(stim(contrast==condition) == stimulus-1  & resp(contrast==condition) == response-1 & conf(contrast==condition)==confidence));
                    end
                end
            end
        end
    end
    
    if compute_measures == 1        
        %Save measures from subject data
        %Across all conditions
        meta_MLE = type2_SDT_MLE(stim,resp,conf,nratings, [], 1);
        mratio_all(subject) = meta_MLE.M_ratio;
        dprime_all(subject) = data_analysis_resp(stim,resp);
        
        %For each individual condition
        for condition = 1:3
            accuracy(subject,condition) = mean(correct(contrast==condition));
            meanConf(subject,condition) = mean(conf(contrast==condition));
            meta_MLE = type2_SDT_MLE(stim(contrast==condition),resp(contrast==condition),conf(contrast==condition),nratings, [], 1);
            mratio(subject,condition) = meta_MLE.M_ratio;
            dprime(subject,condition) = data_analysis_resp(stim(contrast==condition),resp(contrast==condition));
        end
        
    end
    
    %Compute criterion locations
    [~,criteria(subject,:)] = computeSDTcriteria(stim, resp, conf, nratings);
end

%% Save measures from subject data

if save_measures == 1
    exptData.accuracy = accuracy;
    exptData.meanConf = meanConf;
    exptData.mratio = mratio;
    exptData.dprime_all = dprime_all;
    exptData.mratio_all = mratio_all;
    exptData.dprime = dprime;
    exptData.criteria = criteria;
    
    save(['DataForModeling'],'dataCounts','exptData')
end

