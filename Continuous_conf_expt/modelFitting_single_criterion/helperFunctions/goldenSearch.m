function [ metaNoise, fval ] = goldenSearch(a,b,subject)
%Find best fitting metaNoise value
% ------------------------GOLDEN SECTION METHOD----------------------------
% -------------------------------------------------------------------------
% Copyright (c) 2009, Katarzyna Zarnowiec, all rights reserved 
% mailto: katarzyna.zarnowiec@gmail.com
% -------------------------------------------------------------------------

epsilon=0.001;                 % accuracy value
iter= 10;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations

x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);

f_x1=logL_func_metaNoise([x1],subject);       % computing values in x points
f_x2=logL_func_metaNoise([x2],subject);

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1; 
        x1=a+(1-tau)*(b-a); 
        
        f_x1=logL_func_metaNoise([x1],subject);       % computing values in x points
        f_x2=logL_func_metaNoise([x2],subject);
        
    else
        a=x1;
        x1=x2; 
        x2=a+tau*(b-a); 
        
        f_x1=logL_func_metaNoise([x1],subject);       % computing values in x points
        f_x2=logL_func_metaNoise([x2],subject);
        
    end
    
    k=k+1;
end


% chooses minimum point
if(f_x1<f_x2)
    xmin = x1;  fmin = f_x1;
else
    xmin = x2;  fmin = f_x2;
end

metaNoise = xmin; 
fval = fmin; 


end

