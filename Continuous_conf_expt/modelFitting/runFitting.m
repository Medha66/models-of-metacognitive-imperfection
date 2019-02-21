%% Run model fitting 
clear
clc

%Add helper functions
addpath('model_helperFunctions')
global modelToFit
modelToFit = 'normal' %standard, normal, lognormal 

subjects = [1:20];
for subject = subjects
    disp('----------')
    disp(['Fitting subject ' num2str(subject) ', Model: ' modelToFit])
    disp('----------')
    if strcmp(modelToFit,'standard')
        [fval(subject)] = logL_func_standard(subject);
        metaNoise(subject) = 0;
    else
        [metaNoise(subject), fval(subject)] = metaNoise_bestFit(subject);
    end
end

%% Model fits

logL     = -fval;
if strcmp(modelToFit,'standard')
    k    = 34;
else
    k = 35;
end

modelFit.logL = logL;
modelFit.AIC =  2*k - 2*logL 
modelFit.params(:,1) = metaNoise;

save(['fittingResults/fittingResults_' modelToFit], 'modelFit');



