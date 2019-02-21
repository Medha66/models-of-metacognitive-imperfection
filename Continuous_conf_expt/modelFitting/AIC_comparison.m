%% Comparison of AIC values

clear
clc

addpath('fittingResults')

models = {'standard','normal','lognormal'};

for model = 1:3    
    modelName = models{model};
    load(['fittingResults_',modelName])
    
    for subject = 1:20
        AIC(subject,model) = modelFit{subject}.AIC
    end
    
end

% evidence ratio for average difference in AIC values 
evidenceRatio_mean = AICanalysis([mean(AIC)],'e')

% evidence ratio for sum of AIC difference
evidenceRatio_sum = AICanalysis([sum(AIC)],'e')

%% Plot AIC difference b/w lognormal and competing models vs mratio

load DataForModeling.mat
mratio = exptData.mratio_all;
 
subplot(2,2,1)
plot(mratio, AIC(:,3)-AIC(:,1),'o','markersize',7); hold on;
xlabel('meta-d'',d''','fontsize',16)
ylabel('\DeltaAIC')
title('Lognormal vs Standard SDT','fontsize',16)
subplot(2,2,2)
plot(mratio, AIC(:,3)-AIC(:,2),'o','markersize',7);
xlabel('meta-d'',d''','fontsize',16)
ylabel('\DeltaAIC')
title('Lognormal vs Gaussian','fontsize',16)

%AIC difference
AIC_diff(:,1) = AIC(:,3)-AIC(:,1); %Lognormal and standard SDT
AIC_diff(:,2) = AIC(:,3)-AIC(:,2); %Lognormal and Gaussian

% Correlations
[r(1),p(1)] = corr(AIC_diff(:,1),mratio');
[r(2),p(2)] = corr(AIC_diff(:,2),mratio');



