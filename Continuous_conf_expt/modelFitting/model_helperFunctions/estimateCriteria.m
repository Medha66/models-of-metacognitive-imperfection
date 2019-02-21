function [confCriteria, respProb] = estimateCriteria(metaNoise,subject)

% % Estimate best fitting confidence criterion locations for each criterion
% % for given parameters
global modelToFit
global condition

% clear; clc
% global modelToFit
% modelToFit = 'standard';
% metaNoise = 0;
% subject = 4;
% condition = 1;

%Load experiment data
load DataForModeling.mat

%Parameters
%Criteria
c = exptData.criteria(subject,6);
%Criterion locations wrt c
criteria = exptData.criteria(:,:)-exptData.criteria(:,6);
%Remove Inf values
criteria(criteria==Inf) = NaN; criteria(criteria==-Inf) = NaN;

%Express criterion locations as absolute difference from previous criterion
%location
c_negative = abs(wrev(nanmean(criteria(:,1:5)))); c_negative = [c_negative(1), diff(c_negative)];
c_positive = nanmean(criteria(:,7:11)); c_positive = [c_positive(1), diff(c_positive)];
nratings = length(c_positive)+1;

%Lower limit for confidence criteria location parameters (mu)
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
%metaNoise for the Gaussian model (because of crossovers)

if strcmp(modelToFit,'normal')
    w1_ = 1 - normcdf(0,-mu_s1,1);
    w2_ = 1 - normcdf(0,-mu_s2,1);
    
    [q1_]  = evaluateIntegral(-mu_s1,0,metaNoise);
    [q2_]  = evaluateIntegral(-mu_s2,0,metaNoise);
    
    prop_pred_max(1) = (w1_*q1_ + w2_*q2_)/(w1_+w2_);
else
    prop_pred_max(1) = 1;
end

for criterion = 1:5
    %Estimated criterion location on the decision axis
    M = sum(c_negative(1:criterion));
    
    if strcmp(modelToFit,'lognormal')
        %Assuming that the estimated criterion location, M, is the median of the lognormal distribution
        %mu_c1 is the location parameter of the lognormal distribution
        mu_c1(criterion) = log(M);
        
    else %Normal or standard model
        %Criterion location is the mean of the normal distribution
        mu_c1(criterion) = M;
    end
    
    stepSize = .05;
    reversal = 0;
    trial = 1;
    err = 1;
    
    while abs(err) > 10^-4
        
        
        [q1]  = evaluateIntegral(-mu_s1,mu_c1(criterion),metaNoise);
        [q2]  = evaluateIntegral(-mu_s2,mu_c1(criterion),metaNoise);
        
        %Weight q1 and q2 by the total areas that lie to the right of the
        %criterion
        w1 = 1 - normcdf(0,-mu_s1,1);
        w2 = 1 - normcdf(0,-mu_s2,1);
        
        %Predicted proportion of high confidence responses with negative
        %criteria (i.e. for 'left'responses)
        prop_pred(criterion,1) = (w1*q1 + w2*q2)/(w1+w2);
        
        prop_obs(criterion,1) = sum(sum(dataCounts(:,1:6-criterion,subject,condition)))./sum(sum(dataCounts(:,1:nratings,subject,condition)));
        
        
        %% Adjust criterion location
        if prop_obs(criterion,1) < prop_pred_max(1) %Optimization is possible
            err(trial) = prop_obs(criterion,1) - prop_pred(criterion,1);
        else %Optimisation is not possible - set criterion to lowest possible value
            err(trial) = 0; 
            mu_c1(criterion) = lowerLimit; %Lower limit for confidence criteria
        end
        
        if err(trial) > 0 %obs conf > pred conf
            %Decrease confidence criterion
            mu_c1(criterion) = mu_c1(criterion)-stepSize;
            
            if criterion > 1 %Constrain the confidence criteria - lower bound is the previous confidence criterion's location
                if mu_c1(criterion) < mu_c1(criterion-1)
                    mu_c1(criterion) = mu_c1(criterion-1);
                end
            else
                if mu_c1(criterion) < lowerLimit
                    mu_c1(criterion) = lowerLimit;
                    err(trial) = 0;
                end
            end
            direction(trial) = 1;
        elseif err(trial) < 0 %obs conf < pred conf
            %Increase confidence criterion
            mu_c1(criterion) = mu_c1(criterion)+stepSize;
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
    %Save the estimates of - mu of the normal/lognormal distributions
    confCriteria(6-criterion) = mu_c1(criterion);
end

%% Optimise positive confidence criteria

%Calculate the max. prop of high confidence that can be obtained with this
%metaNoise for the Gaussian model

if strcmp(modelToFit,'normal')
    w1_ = 1 - normcdf(0,mu_s1,1);
    w2_ = 1 - normcdf(0,mu_s2,1);
    
    [q1_]  = evaluateIntegral(mu_s1,0,metaNoise);
    [q2_]  = evaluateIntegral(mu_s2,0,metaNoise);
    
    prop_pred_max(2) = (w1_*q1_ + w2_*q2_)/(w1_+w2_);
else
    prop_pred_max(2) = 1;
end

for criterion = 1:5
    M = sum(c_positive(1:criterion));
    
    if strcmp(modelToFit,'lognormal')
        %Assuming that the criterion location is median of the lognormal distribution
        mu_c2(criterion) = log(M);
    else
        %Criterion location is the mean of the normal distribution
        mu_c2(criterion) = M;
    end
    
    stepSize = .05;
    reversal = 0;
    trial = 1;
    err = 1;
    
    while abs(err) > 10^-4
        %Assuming that the criterion location is median of the distribution
        
        [q1] = evaluateIntegral(mu_s1,mu_c2(criterion),metaNoise);
        [q2] = evaluateIntegral(mu_s2,mu_c2(criterion),metaNoise);
        
        %Weight q1 and q2 by the total areas that lie to the right of the
        %criterion
        w1 = 1 - normcdf(0,mu_s1,1);
        w2 = 1 - normcdf(0,mu_s2,1);
        
        %Predicted proportion of high confidence responses with positive
        %criteria (i.e. for 'right'responses)
        prop_pred(criterion,2) = (w1*q1 + w2*q2)/(w1+w2);
        
        prop_obs(criterion,2) = sum(sum(dataCounts(:,6+criterion+1:12,subject,condition)))./sum(sum(dataCounts(:,nratings+1:12,subject,condition)));
        
        
        %% Adjust criterion location
        if prop_obs(criterion,2) < prop_pred_max(2)
            err(trial) = prop_obs(criterion,2) - prop_pred(criterion,2);
        else
            err(trial) = 0;
            mu_c2(criterion) = lowerLimit;
        end
        
        if err(trial) > 0 %obs conf > pred conf
            %Decrease confidence criterion
            mu_c2(criterion) = mu_c2(criterion)-stepSize;
            direction(trial) = 1;
            
            if criterion > 1
                if mu_c2(criterion) < mu_c2(criterion-1)
                    mu_c2(criterion) = mu_c2(criterion-1);
                end
            else
                if mu_c2(criterion) < lowerLimit
                    mu_c2(criterion) = lowerLimit;
                    err(trial) = 0;
                end
            end
            
        elseif err(trial) < 0  %obs conf < pred conf
            %Increase confidence criterion
            mu_c2(criterion) = mu_c2(criterion)+stepSize;
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
    
    confCriteria(6+criterion) = mu_c2(criterion);
end


%% Get response probabilities (prob of high conf) from the criteria

for criterion = 1:5
    %Left confidence criteria
    [~,respProb_(1,6-criterion),total_area(1,1)] = evaluateIntegral(-mu_s1,mu_c1(criterion),metaNoise);
    [~,respProb_(2,6-criterion),total_area(2,1)] = evaluateIntegral(-mu_s2,mu_c1(criterion),metaNoise);
    
    %Right confidence criteria
    [~,respProb_(1,5+criterion),total_area(1,2)] = evaluateIntegral(mu_s1,mu_c2(criterion),metaNoise);
    [~,respProb_(2,5+criterion),total_area(2,2)] = evaluateIntegral(mu_s2,mu_c2(criterion),metaNoise);
end

%% cdf of the stimulus distributions at 0

for stim = 1:2
    respProb(stim,1:6) = diff([0 respProb_(stim,1:5) total_area(stim,1)])./2; %Left confCriteria
    respProb(stim,7:12) = -diff([total_area(stim,2) respProb_(stim,6:10) 0])./2; %Right confCriteria
end



