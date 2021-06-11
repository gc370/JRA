function [Eulersmethodresults,Eulersimprovedmethodresults,Midpointmethodresults] = callDE(stepsize,xinitial,tinitial,tend,xprime)

% function handles for inserting equations 

%                            x' =  @(x) x^2-6 e.g
%                            x' = @(x,t) = x^2 - t e.g


[tavlues,xvaluesmidpoint] = Midpoint(stepsize,xinitial,tinitial,tend,xprime);

[tvalues,xvalueseuler] = ImprovedEuler(stepsize,xinitial,tinitial,tend,xprime);

[tvalues,xvalues] = Eulersmethod(stepsize,xinitial,tinitial,tend,xprime);

[tvalues,xrungekuttavalues] = RungeKuttaOrder4(stepsize,xinitial,tinitial,tend,xprime);


%Eulersmethodresults = transpose([tvalues;xvalues]);
%Eulersimprovedmethodresults = transpose([tvalues;xvalueseuler]);
%Midpointmethodresults = transpose([tvalues;xvaluesmidpoint]);

mergedatafortable = transpose([tvalues; xvalues; xvalueseuler; xvaluesmidpoint;xrungekuttavalues]);
table = array2table(mergedatafortable, 'VariableNames',{'Time t','Eulers Method','Improved Eulers Method','Midpoint Method','Runge Kutta Order 4'})




end

