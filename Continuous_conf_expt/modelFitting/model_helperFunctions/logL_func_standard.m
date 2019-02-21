function logL = logL_func_standard(subject)
%Computes the log likelihood associated with the standard model

load DataForModeling
global condition

logL = 0;

for condition = 1:3
    [confCriteria, respProb] = estimateCriteria(0,subject);
    for stim = 0:1                
        respProb(respProb==0)=10^-10;
        logL = logL - sum(log(respProb(stim+1,:)) .* dataCounts(stim+1,:,subject,condition));
        
    end
end
end

