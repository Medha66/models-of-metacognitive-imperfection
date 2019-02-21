function [metaNoise, fval] = metaNoise_bestFit(subject)

%Function to return logL and best fitting metaNoise 

% Initial fitting - plot logL for different values of metaNoise

% global criterion
%global condition
% subject = 4;
% criterion = 5;
% condition = 1;

global modelToFit

clear metaNoise logL

if strcmp(modelToFit,'lognormal')
    lowerLimit = .1;
else
    lowerLimit = .05;
end

%Initial estimate
metaNoise_vect = [lowerLimit:.2:3];

for ind = 1:length(metaNoise_vect)    
    logL(ind) = logL_func_metaNoise([metaNoise_vect(ind)],subject);
end

% % Plot logL vs metaNoise
% figure; hold on; 
% plot(metaNoise_vect,logL,'o-k')
% title(['subject:',num2str(subject)])

%Find min logL values
[~,index] = sort(logL);
metaNoise_init = metaNoise_vect(index(1));
 
%Final estimate - golden search
a = metaNoise_init - .15; %Lower bound
b = metaNoise_init + .15; %Upper bound
a(a<lowerLimit)=lowerLimit;

for iter = 1:1
    [metaNoise_(iter),fval_(iter)] = goldenSearch(a(iter),b(iter),subject);
end

[fval,minIndex] = min(fval_);
metaNoise = metaNoise_(minIndex);



end
