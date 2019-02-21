function [metaNoise,fval] = metaNoise_bestFit(subject)

% Function to return logL and best fitting metaNoise 

% Initial fitting - plot logL for different values of metaNoise

% clear 
% clc
% subject = 4;
% condition = 1;

global modelToFit

if strcmp(modelToFit,'lognormal')
    lowerLimit = .05;
else
    lowerLimit = .05;
end

%Coarse search
metaNoise_vect = [lowerLimit:.1:2.5];
for ind = 1:length(metaNoise_vect)
    logL(ind) = logL_func_metaNoise([metaNoise_vect(ind)],subject);
end


[~,minInd] = min(logL);

%Plot
%figure;
plot(metaNoise_vect,logL,'o-','color',[rand rand rand],'linewidth',1.5); hold on
xlabel metaNoise
ylabel logL
title(['subject: ',num2str(subject)])

%Find best coasrse estimate of metaNoise from the logL function 
 metaNoise_init = metaNoise_vect(minInd);

 
%Final estimate - golden search
a = metaNoise_init - .15; %Lower bound
b = metaNoise_init + .15; %Upper bound
a(a<lowerLimit)=lowerLimit;

[metaNoise,fval] = goldenSearch(a,b,subject);


end
