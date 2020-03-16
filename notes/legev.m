%function y=legev(x,range,order);
% x:scalr
% range: vector two elements
% order: Degree for the polynomial evaluation (scalar)
%
% y: scalar
%
%
%
%
% Chebysheb polinomial of degree (order) evaluated at point (x), where the range scales to map in
% closed interval [-1 1]

function y=legev(x,range,order);

if isnan(x)==1
   y=NaN;
else
   
   
   
   %Normalization for the input range to [-1,1]
   xmax=max(range);xmin=min(range);
   x1=2/(xmax-xmin)*x-(xmax+xmin)/(xmax-xmin);

   y=cos(order*acos(x1));
   
   
   %y=pp;
   %y=sqrt((2*order+1)/2)*pp;
   
end


