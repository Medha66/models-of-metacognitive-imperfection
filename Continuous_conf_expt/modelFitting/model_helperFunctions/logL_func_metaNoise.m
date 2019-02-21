function logL = logL_func_metaNoise(metaNoise,subject)
%Computes the log likelihood associated with the specified value of metaNoise

load DataForModeling
global condition

logL = 0;

for condition = 1:3
    %Find the best fitting confidence criterion locations for given
    %metaNoise
    [confCriteria, respProb] = estimateCriteria(metaNoise,subject);
    
    for stim = 0:1
        
        respProb(respProb==0)=10^-10;
        logL = logL - sum(log(respProb(stim+1,:)) .* dataCounts(stim+1,:,subject,condition));
        
    end
end




