function [tvalues, xvaluesmidpoint, xprimevalues] = Midpoint(stepsize,xinitial,tinitial,tend,xprime)

% --------- Midpoint Method for approximating ODE's 

% function handles for inserting equations 

%                            x' =  @(x) x^2-6 e.g
%                            x' = @(x,t) = x^2 - t e.g

% --------- General derivation 

% x(n+1) = x(n) + hf(x(n) + (h/2)*f(X(n),t(n)))
% T values are just counting the step size incremements
% Xprime values are calulcated from subbing in x and t into x'.

if(tinitial > tend)
   disp("Please make sure we are going forward in time :)");
else

%Making our DE into a discretized/recursion format ready for Eulers    
    
xprime = string(char(xprime));
xprime = replace(xprime,["x","t"],["xvaluesmidpoint(counter)","tvalues(counter)"]);
xprime = replace(xprime,"@(xvaluesmidpoint(counter),tvalues(counter))","@(xvaluesmidpoint,tvalues) ");
xprime = replace(xprime,"@(tvalues(counter),xvalues(counter))","@(xvaluesmidpoint,tvalues)");
xprime = replace(xprime,"@(xvaluesmidpoint,tvalues)","@(xvaluesmidpoint,tvalues,counter)");
xprime = str2func(xprime);

%input Values

h = stepsize;
xvaluesmidpoint = [];
xprimevalues = [];
tvalues = [];
counter = 1;
maxcounter = floor((tend-tinitial)/h)+1;
tvalues(1) = tinitial;
xvaluesmidpoint(1) = xinitial;
xprimevalues(counter) = xprime(xvaluesmidpoint(counter),tvalues(counter),1);

while (counter < maxcounter && counter < 10000)
  
   tvalues(counter+1) = tvalues(1) + h*counter;
   xvaluesmidpoint(counter+1) = xvaluesmidpoint(counter) + h*(xprime(xvaluesmidpoint(counter) + (h/2)*xprime(xvaluesmidpoint(counter),tvalues(counter),1),tvalues(counter) + (h/2),1));
   xprimevalues(counter+1) = xprime(xvaluesmidpoint(counter+1),tvalues(counter+1),1);
   counter = counter + 1;
   
end

%plot(tvalues,xvaluesmidpoint)

end