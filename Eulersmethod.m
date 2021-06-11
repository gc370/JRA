function [tvalues,xvalues,xprimevalues] = Eulersmethod(stepsize,xinitial,tinitial,tend,xprime)

% --------- Eulers Method for approximating ODE's 

% function handles for inserting equations 

%                            x' =  @(x) x^2-6 e.g
%                            x' = @(x,t) = x^2 - t e.g

% --------- General derivation 

%   x(t + h) = x(t) + hf (t, x(t)) + R1 (t) -> x(t) + hx'(t) as h -> sufficiently small.

if(tinitial > tend)
   disp("Please make sure we are going forward in time :)");
else

%Making our DE into a discretized/recursion format ready for Eulers    
    
xprime = string(char(xprime));
xprime = replace(xprime,["x","t"],["xvalues(counter)","tvalues(counter)"]);
xprime = replace(xprime,"@(xvalues(counter),tvalues(counter))","@(xvalues,tvalues) ");
xprime = replace(xprime,"@(tvalues(counter),xvalues(counter))","@(xvalues,tvalues)");
xprime = replace(xprime,"@(xvalues,tvalues)","@(xvalues,tvalues,counter)");
xprime = str2func(xprime);

%input Values

h = stepsize;
xvalues = [];
xprimevalues = [];
tvalues = [];
counter = 1;
maxcounter = floor((tend-tinitial)/h)+1;
tvalues(1) = tinitial;
xvalues(1) = xinitial;
xprimevalues(counter) = xprime(xvalues(counter),tvalues(counter),1);

% xprimevalues(counter) = xprime(xvalues(counter),tvalues(counter),1); need to find out 
% why 1 is needed for counter)

while (counter < maxcounter && counter < 10000)
  
   tvalues(counter+1) = tvalues(1) + h*counter;
   xvalues(counter+1) = xvalues(counter) + h*(xprimevalues(counter));
   xprimevalues(counter+1) = xprime(xvalues(counter+1),tvalues(counter+1),1);
   counter = counter + 1;
   
end

%plot(tvalues,xvalues)

end