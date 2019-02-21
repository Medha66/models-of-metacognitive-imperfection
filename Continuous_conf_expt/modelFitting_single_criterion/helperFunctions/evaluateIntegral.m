function [q, q1,q2] = evaluateIntegral(mu1,mu2,metaNoise)
% Numerically evaluate the double integral with lognormal and normal functions
% to determine the prop of high confidence trials
%Location parameter - mu2, scale parameter - metaNoise

sigma = 1;
global modelToFit

if strcmp(modelToFit,'lognormal')
    %Lognormal
    
    %Area of high confidence
    fun1 = @(y,x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2./2)).*(1./(sqrt(2.*pi).*y*metaNoise).*exp(-(log(y)-mu2).^2./(2.*metaNoise.^2)));
    
    ymin = 0;
    ymax = Inf;
    xmin = @(y) y;
    xmax = Inf;
    q1 = integral2(fun1,ymin,ymax,xmin,xmax);
    %q1 = integral2(fun1,ymin,ymax,xmin,xmax,'AbsTol',10e-50,'RelTol',10e-10);
    
    %Different way of computing the same integral
    % xmin = 0;
    % xmax = Inf;
    % ymin = 0;
    % ymax = @(x) x;
    % q1 = integral2(fun1,xmin,xmax,ymin,ymax);
    
    %Total area
    fun2 = @(x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2/(2.*sigma^2)));
    
    xmin = 0;
    xmax = Inf;
    q2 = integral(fun2,xmin,xmax);
    
    q = q1/q2;
    
else
    %Gaussian additive noise    
    fun1 = @(x,y) (1./sqrt(2.*pi).*exp(-(x-mu1).^2./2)).*(1./(sqrt(2.*pi).*metaNoise).*exp(-(y-mu2).^2./(2.*metaNoise.^2)));
    
    %High conf
    ymin = -Inf;
    ymax = @(x) x;
    xmin = 0;
    xmax = Inf;
    q1 = integral2(fun1,xmin,xmax,ymin,ymax);    
    
    %Total area
    fun2 = @(x) (1./sqrt(2.*pi).*exp(-(x-mu1).^2/(2.*sigma^2)));
    
    xmin = 0;
    xmax = Inf;
    q2 = integral(fun2,xmin,xmax);
    
    q = q1/q2;           
    
end



end