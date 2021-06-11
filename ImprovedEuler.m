function [tvalues, xvalueseuler, xprimevalues] = ImprovedEuler(stepsize,xinitial,tinitial,tend,xprime)

% --------- Midpoint Method for approximating ODE's 

% function handles for inserting equations 

%                            x' =  @(x) x^2-6 e.g
%                            x' = @(x,t) = x^2 - t e.g

% --------- General derivation 

% We take the average gradient between to points (t_n,x_n) and
% (t_n+1,x_n+1) where x_n+1 is approximate using Eulers scheme.

% I.e we have ******X_n+1 = X_n + h*x'(n)****** (EULERS)

% Then using the midpoint method we have -> 

% X_n+1 = X_n + 0.5*h*(x'(n) + x'(n+1)) where x'(n+1) uses Eulers approx.

if(tinitial > tend)
   disp("Please make sure we are going forward in time :)");
else

%Making our DE into a discretized/recursion format ready for Eulers    
    
xprime = string(char(xprime));
xprime = replace(xprime,["x","t"],["xvalueseuler(counter)","tvalues(counter)"]);
xprime = replace(xprime,"@(xvalueseuler(counter),tvalues(counter))","@(xvalueseuler,tvalues) ");
xprime = replace(xprime,"@(tvalues(counter),xvalueseuler(counter))","@(xvalueseuler,tvalues)");
xprime = replace(xprime,"@(xvalueseuler,tvalues)","@(xvalueseuler,tvalues,counter)");
xprime = str2func(xprime);

%input Values

h = stepsize;
xvalueseuler = [];
xprimevalues = [];
tvalues = [];
counter = 1;
maxcounter = floor((tend-tinitial)/h)+1;
tvalues(1) = tinitial;
xvalueseuler(1) = xinitial;
xprimevalues(counter) = xprime(xvalueseuler(counter),tvalues(counter),1);



while (counter < maxcounter && counter < 10000)
  
   tvalues(counter+1) = tvalues(1) + h*counter;
   xvalueseuler(counter+1) = xvalueseuler(counter) + (h/2)*(xprimevalues(counter) + xprime(xvalueseuler(counter)+ h*xprimevalues(counter),tvalues(counter+1),1));
   xprimevalues(counter+1) = xprime(xvalueseuler(counter+1),tvalues(counter+1),1);
   counter = counter + 1;
   
end

%plot(tvalues,xvalueseuler);

end