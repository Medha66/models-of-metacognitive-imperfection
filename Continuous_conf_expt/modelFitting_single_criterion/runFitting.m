%% Run fitting for all subjects
%This model fitting procedure estimates metaNoise separately for each
%criterion location
clear
clc

global modelToFit
global condition
global criterion


modelToFit = 'lognormal' %normal or lognormal

subjects = [1:20];

for criterion = 1:5
    for condition = 1:3
        for subject = subjects
            disp('----------')
            disp(['Fitting subject ' num2str(subject) '     Contrast ' num2str(condition),'   Criterion ' num2str(criterion),'    Model: ' modelToFit])
            disp('----------')
            
            [metaNoise(subject,condition,criterion), fval(subject,condition,criterion)] = metaNoise_bestFit(subject);
            
        end
    end
end

%% Model fits

logL     = -fval;
if strcmp(modelToFit,'standard')
    k    = 30;
else
    k    = 31;
end

modelFit.logL = logL;
modelFit.AIC = -2*logL + 2*k;
modelFit.params = metaNoise;

save(['fittingResults/fittingResults_' modelToFit], 'modelFit');


