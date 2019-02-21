function [ X ] = correlatedRV(p1, p2, model, N)
%Generate perfectly correlated random variables from
%multivariate normal random variables (copula method)

%p1 is the vector containing mean of each criterion 
%p2 is metaNoise parameter

%Generate multivariate normal RVs of mean 0 and SD 1
%each distribution corresponding to one confidence criterion

ncriteria = length(p1);
CovMat = ones([ncriteria ncriteria]); %Perfect correlation between the RVs
Z = mvnrnd(zeros([1,length(p1)]), CovMat, N);
U = normcdf(Z);

if strcmp(model,'normal')
    for crit = 1:ncriteria
        X(crit,:) = norminv(U(:,crit),p1(crit),p2);
    end
elseif strcmp(model,'lognormal')
    for crit = 1:ncriteria
        X(crit,:) = logninv(U(:,crit),p1(crit),p2);
    end
elseif strcmp(model,'gamma')
    for crit = 1:ncriteria
        X(crit,:) = gaminv(U(:,crit),p1(crit),p2);
    end
end

