% function [confCriteria, respProb] = estimateCriteria(metaNoise,subject)
% 
% % % Estimate best fitting confidence criterion locations for each criterion
% % location
% global modelToFit
% global condition
% global criterion %From 1 to 5

clear; clc
global modelToFit
metaNoise = .9;
subject = 1;
condition = 1;
criterion = 1;
modelToFit = 'normal';

%Load experiment data
load DataForModeling.mat

%Parameters
%Criteria
c = exptData.criteria(subject,6);
criteria = exptData.criteria(:,:)-exptData.criteria(:,6);
%Remove Inf values
criteria(criteria==Inf) = NaN; criteria(criteria==-Inf) = NaN;

%Criterion locations wrt c
c_negative = abs(wrev(nanmean(criteria(:,1:5)))); c_negative = [c_negative(1), diff(c_negative)];
c_positive = nanmean(criteria(:,7:11)); c_positive = [c_positive(1), diff(c_positive)];
nratings = length(c_positive)+1;

%Lower limit for confidence criteria
if strcmp(modelToFit,'lognormal')
    lowerLimit = -6;
else
    lowerLimit = 10^-5;
end

%Stimulus distributions for s1 and s2 - shift the mean of the two
%distributions such that the origin is defined at c
dprime = exptData.dprime(subject,condition);
mu_s1 =  -dprime/2-c;
mu_s2 =  dprime/2-c;

%% Optimise negative confidence criteria

%Calculate the max. prop of high confidence that can be obtained with this
%metaNoise

M = sum(c_negative(1:criterion));

if strcmp(modelToFit,'normal')
    w1_ = 1 - normcdf(0,-mu_s1,1);
    w2_ = 1 - normcdf(0,-mu_s2,1);
    
    [q1_]  = evaluateIntegral(-mu_s1,0,metaNoise);
    [q2_]  = evaluateIntegral(-mu_s2,0,metaNoise);
    
    %Criterion location is the mean of the normal distribution
    mu_c1 = M;
    prop_pred_max(1) = (w1_*q1_ + w2_*q2_)/(w1_+w2_);
    
else 
    prop_pred_max(1) = 1;
    %Assuming that the criterion location is median of the lognormal distribution
    mu_c1 = log(M);
end

stepSize = .05;
reversal = 0;
trial = 1;
err = 1;

while abs(err) > 10^-4
    
    %Proportion of high confidence normalized by proportion of responses 
    %From S1 distribution (hits)
    [q1]  = evaluateIntegral(-mu_s1,mu_c1,metaNoise);
    %From S2 distribution (false alarms)
    [q2]  = evaluateIntegral(-mu_s2,mu_c1,metaNoise);
    
    %Weight q1 and q2 by the total areas that lie to the right of the
    %criterion
    w1 = 1 - normcdf(0,-mu_s1,1);
    w2 = 1 - normcdf(0,-mu_s2,1);
    
    %Predicted proportion of high confidence responses with negative
    %criteria (i.e. for 'left'responses)
    prop_pred(1) = (w1*q1 + w2*q2)/(w1+w2);
    
    prop_obs(1) = sum(sum(dataCounts(:,1:6-criterion,subject,condition)))./sum(sum(dataCounts(:,1:nratings,subject,condition)));
    
    
    %% Adjust criterion location
    if prop_obs(1) < prop_pred_max(1)
        err(trial) = prop_obs(1) - prop_pred(1);
    else
        err(trial) = 0;
        mu_c1 = lowerLimit; %Lower limit for confidence criteria
    end
    
    if err(trial) > 0 %obs conf > pred conf
        %Decrease confidence criterion
        mu_c1 = mu_c1-stepSize;
        
        
        if mu_c1 < lowerLimit
            mu_c1 = lowerLimit;
            err(trial) = 0;
        end
        
        direction(trial) = 1;
    elseif err(trial) < 0 %obs conf < pred conf
        %Increase confidence criterion
        mu_c1 = mu_c1+stepSize;
        direction(trial) = -1;
    end
    
    %Check if loop is stuck
    if trial > 1 && err(trial) == err(trial-1)
        display breaking
        break
    end
    
    %Decrease stepsize at reversal
    if trial > 1 && direction(trial) ~= direction(trial-1)
        reversal = reversal + 1;
        stepSize = stepSize/2;
    end
    trial = trial+1;
end

confCriteria(1) = mu_c1;

%% Optimise positive confidence criteria

M = sum(c_positive(1:criterion));

if strcmp(modelToFit,'normal')
    w1_ = 1 - normcdf(0,mu_s1,1);
    w2_ = 1 - normcdf(0,mu_s2,1);
    
    [q1_]  = evaluateIntegral(mu_s1,0,metaNoise);
    [q2_]  = evaluateIntegral(mu_s2,0,metaNoise);
    
    prop_pred_max(2) = (w1_*q1_ + w2_*q2_)/(w1_+w2_);
    
    %Criterion location is the mean of the normal distribution
    mu_c2 = M;
    
elseif strcmp(modelToFit,'lognormal')
    prop_pred_max(2) = 1;
    %Assuming that the criterion location is median of the lognormal distribution
    mu_c2 = log(M);
end

stepSize = .05;
reversal = 0;
trial = 1;
err = 1;

while abs(err) > 10^-4
    %Assuming that the criterion location is median of the distribution
    
    %Proportion of high confidence normalized by proportion of responses 
    %From S1 distribution (false alarms)
    [q1] = evaluateIntegral(mu_s1,mu_c2,metaNoise);
    %From S2 distribution (hits)
    [q2] = evaluateIntegral(mu_s2,mu_c2,metaNoise);
    
    %Weight q1 and q2 by the total areas that lie to the right of the
    %criterion
    w1 = 1 - normcdf(0,mu_s1,1);
    w2 = 1 - normcdf(0,mu_s2,1);
    
    %Predicted proportion of high confidence responses with positive
    %criteria (i.e. for 'right'responses)
    prop_pred(2) = (w1*q1 + w2*q2)/(w1+w2);
    
    prop_obs(2) = sum(sum(dataCounts(:,6+criterion+1:12,subject,condition)))./sum(sum(dataCounts(:,nratings+1:12,subject,condition)));
    
    
    %% Adjust criterion location
    if prop_obs(2) < prop_pred_max(2)
        err(trial) = prop_obs(2) - prop_pred(2);
    else
        err(trial) = 0;
        mu_c2 = lowerLimit;
    end
    
    if err(trial) > 0 %obs conf > pred conf
        %Decrease confidence criterion
        mu_c2 = mu_c2-stepSize;
        direction(trial) = 1;
        
        
        if mu_c2 < lowerLimit
            mu_c2 = lowerLimit;
            err(trial) = 0;
        end
        
        
    elseif err(trial) < 0  %obs conf < pred conf
        %Increase confidence criterion
        mu_c2 = mu_c2+stepSize;
        direction(trial) = -1;
    end
    
    %Check if loop is stuck
    if trial > 1 && err(trial) == err(trial-1)
        display breaking
        break
    end
    
    %Decrease stepsize at reversal
    if trial > 1 && direction(trial) ~= direction(trial-1)
        reversal = reversal + 1;
        stepSize = stepSize/2;
    end
    
    trial = trial+1;
    
end

confCriteria(3) = mu_c2;


%% Get response probabilities (prob of high conf) from the criteria

%Left confidence criteria
[~,respProb_(1,1),total_area(1,1)] = evaluateIntegral(-mu_s1,mu_c1,metaNoise);
[~,respProb_(2,1),total_area(2,1)] = evaluateIntegral(-mu_s2,mu_c1,metaNoise);

%Right confidence criteria
[~,respProb_(1,2),total_area(1,2)] = evaluateIntegral(mu_s1,mu_c2,metaNoise);
[~,respProb_(2,2),total_area(2,2)] = evaluateIntegral(mu_s2,mu_c2,metaNoise);

%% cdf of the stimulus distributions at 0

for stim = 1:2
    respProb(stim,1:2) = diff([0 respProb_(stim,1) total_area(stim,1)])./2; %Left confCriteria
    respProb(stim,3:4) = -diff([total_area(stim,2) respProb_(stim,2) 0])./2; %Right confCriteria
    
    %Just to check if resp probabilities match the data
    counts(stim,:) = [sum(dataCounts(stim,1:nratings-criterion,subject,condition)),sum(dataCounts(stim,nratings+1-criterion:nratings,subject,condition))...
        ,sum(dataCounts(stim,nratings+1:nratings+criterion,subject,condition)),sum(dataCounts(stim,nratings+1+criterion:end,subject,condition))];
end



