function [Approximationsofu] = Finitedifferencessecondorder(tinitial,tend,xtinitial,xtend,xdoubleprime,Noofpoints)


points = Noofpoints;
h = ((tend-tinitial)/points)
pointx = [];
for i = 1: points+1;
pointx(i) = tinitial + (i-1)*h;
end

% h = D(x) Delta X ------ Now we make our Interval of points.
if (points < 3)
    Disp("Please add more points)")
else


% Setting up RHS for (Delta x)^2f(xi) where f is x''.

RHS = zeros(1,points+1);
RHS(1) = xtinitial;
for i = 2: points
RHS(i) = (h^2)*xdoubleprime(pointx(i));
end
RHS(end) = xtend;
RHS = transpose(RHS);

%setting up the left hand side matrix.

LHS = zeros(points+1,points+1);
LHS(1,1)=1;
LHS(points+1,points+1)=1;
for i = 2: points
    LHS(i,i) = -2;
    LHS(i,i-1) = 1;
    LHS(i,i+1) =1;
end


% We now show the approximate values of the function U at the intervals
% X_i + (i-1)*h

LHS = inv(LHS);
Approximationsofu = LHS * RHS;

for i = 1: points+1
    disp("The approximate value at x = " + pointx(i) + " is " + Approximationsofu(i))
end

end

