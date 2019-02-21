%% Plot logL against metaNoise

subject = 2;
global condition
global modelToFit 

modelToFit = 'lognormal'

for condition = 1:3

metaNoise_vect = [.5:.2:2];

for ind = 1:length(metaNoise_vect)
    
    logL_init(ind) = logL_func_metaNoise([metaNoise_vect(ind)],subject);
end

% Plot logL vs metaNoise
figure; hold on; 
plot(metaNoise_vect,logL_init,'o-k')
title(['subject:',num2str(subject)])
end

