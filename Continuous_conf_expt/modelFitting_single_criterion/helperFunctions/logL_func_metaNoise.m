function logL = logL_func_metaNoise(metaNoise,subject)
%Estimate the log likelihood associated with the parameters of the model

load DataForModeling
global condition
global criterion

[confCriteria, respProb] = estimateCriteria2(metaNoise,subject);

nratings = length(dataCounts(1,:,1,1))/2;

for stim = 0:1
    counts(stim+1,:) = [sum(dataCounts(stim+1,1:nratings-criterion,subject,condition)),sum(dataCounts(stim+1,nratings+1-criterion:nratings,subject,condition))...
        ,sum(dataCounts(stim+1,nratings+1:nratings+criterion,subject,condition)),sum(dataCounts(stim+1,nratings+1+criterion:end,subject,condition))];
end

logL = 0;

respProb(respProb==0)=10^-10;
logL = logL - sum(sum((log(respProb) .* counts)));


end

