function [tvalues,xrungekuttavalues,xprimevalues] = RungeKuttaOrder4(stepsize,xinitial,tinitial,tend,xprime)
% function handles for inserting equations 

%                            x' =  @(x) x^2-6 e.g
%                            x' = @(x,t) = x^2 - t e.g


% ------ RUNGE KUTTA ORDER 4 -------

% --- Formula

%   k1 = hf(x(i),t(i))

%   k2 = hf(x(i)+ (0.5*k1),t(i)+ h/2)

%   k3 = hf(x(i) +(0.5*k2),t(i)+ h/2)

%   k4 = hf(x(i) + k3, t(i+1)

%   X(i+1) = x(i) + 1/6(k1 + 2k2 + 2k3 + k4)


if(tinitial > tend)
   disp("Please make sure we are going forward in time :)");
else

%Making our DE into a discretized/recursion format ready for Eulers    
    
xprime = string(char(xprime));
xprime = replace(xprime,["x","t"],["xrungekuttavalues(counter)","tvalues(counter)"]);
xprime = replace(xprime,"@(xrungekuttavalues(counter),tvalues(counter))","@(xrungekuttavalues,tvalues) ");
xprime = replace(xprime,"@(tvalues(counter),xvalues(counter))","@(xrungekuttavalues,tvalues)");
xprime = replace(xprime,"@(xrungekuttavalues,tvalues)","@(xrungekuttavalues,tvalues,counter)");
xprime = str2func(xprime);

%input Values

k1 = [];
k2 = [];
k3 = [];
k4 = [];

h = stepsize;
xrungekuttavalues = [];
xprimevalues = [];
tvalues = [];
counter = 1;
maxcounter = floor((tend-tinitial)/h)+1;
tvalues(1) = tinitial;
xrungekuttavalues(1) = xinitial;
xprimevalues(counter) = xprime(xrungekuttavalues(counter),tvalues(counter),1);

while (counter < maxcounter && counter < 10000)
  
   tvalues(counter+1) = tvalues(1) + h*counter;
   k1(counter) = h*xprime(xrungekuttavalues(counter),tvalues(counter),1);
   k2(counter) = h*xprime(xrungekuttavalues(counter) + (0.5*k1(counter)),tvalues(counter) + (h/2),1);
   k3(counter) = h*xprime(xrungekuttavalues(counter) + (0.5*k2(counter)),tvalues(counter) + (h/2),1);
   k4(counter) = h*xprime(xrungekuttavalues(counter) + k3(counter),tvalues(counter+1),1);
   xrungekuttavalues(counter+1) = xrungekuttavalues(counter) + (1/6)* (k1(counter) + 2*k2(counter) + 2*k3(counter) + k4(counter));
   xprimevalues(counter+1) = xprime(xrungekuttavalues(counter+1),tvalues(counter+1),1);
   counter = counter + 1;
   
end

%plot(tvalues, xrungekuttavalues);

end


  
