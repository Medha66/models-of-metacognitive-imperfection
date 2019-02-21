function [q, q1,q2] = evaluateIntegral(mu1,mu2,metaNoise)
% Numerically evaluate the double integral with gamma and normal confidence
% criterion distributions and Gaussian evidence distributions
% to determine the prop of high confidence trials
%Conf criterion distribution - location parameter - mu2, scale parameter - metaNoise
%Evidence distribution - mu1 - mean

sigma = 1; %SD of evidence distributions
global modelToFit

%Area of high conf
if strcmp(modelToFit,'lognormal')
    %Lognormal
    
    %Area of high confidence
    fun1 = @(y,x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2./2)).*(1./(sqrt(2.*pi).*y*metaNoise).*exp(-(log(y)-mu2).^2./(2.*metaNoise.^2)));
    
    ymin = 0;
    ymax = Inf;
    xmin = @(y) y;
    xmax = Inf;
    q1 = integral2(fun1,ymin,ymax,xmin,xmax);
    %For more precision -
    %q1 = integral2(fun1,ymin,ymax,xmin,xmax,'AbsTol',10e-50,'RelTol',10e-10);
    
    %Alternate way of computing the integral
    % xmin = 0;
    % xmax = Inf;
    % ymin = 0;
    % ymax = @(x) x;
    % q1 = integral2(fun1,xmin,xmax,ymin,ymax);
    
    %Type-1 response frequency associated with the distribution
    fun2 = @(x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2/(2.*sigma^2)));
    
    xmin = 0;
    xmax = Inf;
    q2 = integral(fun2,xmin,xmax);
    
    %Normalized area of "high confidence"
    q = q1/q2;
    
elseif strcmp(modelToFit,'normal')
    %Gaussian additive noise
    fun1 = @(x,y) (1./sqrt(2.*pi).*exp(-(x-mu1).^2./2)).*(1./(sqrt(2.*pi).*metaNoise).*exp(-(y-mu2).^2./(2.*metaNoise.^2)));
    %fun1 = @(y,x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2./2)).*(1./(sqrt(2.*pi).*metaNoise).*exp(-(y-mu2).^2./(2.*metaNoise.^2)));
    
        %High conf
        ymin = -Inf;
        ymax = @(x) x;
        xmin = 0;
        xmax = Inf;
        q1 = integral2(fun1,xmin,xmax,ymin,ymax);
    
    %Alternate way of computing the integral
%     ymin = -Inf;
%     ymax = Inf;
%     xmin = @(y) y;
%     xmax = Inf;
%     q1 = integral2(fun1,ymin,ymax,xmin,xmax);
    
    %Total area
    fun2 = @(x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2/(2.*sigma^2)));
    
    xmin = 0;
    xmax = Inf;
    q2 = integral(fun2,xmin,xmax);
    
    q = q1/q2;
    
elseif strcmp(modelToFit,'standard') %Standard model with no metaNoise
    %Area of high confidence
    fun2 = @(x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2/(2.*sigma^2)));
    
    xmin = mu2;
    xmax = Inf;
    q1 = integral(fun2,xmin,xmax);
    
    %Total area
    xmin = 0;
    xmax = Inf;
    q2 = integral(fun2,xmin,xmax);
    
    q = q1/q2;
    
    
end



end