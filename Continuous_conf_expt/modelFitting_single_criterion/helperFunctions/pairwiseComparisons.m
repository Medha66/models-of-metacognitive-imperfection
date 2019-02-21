function [tstat,meanDiff,pval ] = pairwiseComparisons(y)
%Conduct paired t-tests - criterion 1 against others and get with sbject
%sem

%t-tests
for test = 2:5
    [~,pval(test-1),~,stats] = ttest(y(:,1),y(:,test));
    tstat(test-1,1) = stats.tstat;
    meanDiff(test-1,1) = mean(y(:,1))-mean(y(:,test));
end

end

