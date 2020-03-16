%function Y=legev_vec(X,range,order)
% X: Vector of elements
% range: Vector 2 elements
% order: scalar
%
% Y: vector
%
%
%
% Chebysheb polinomial of degree (order) evaluated at points (X), where the range scales to map in
% closed interval [-1 1]


function Y=legev_vec(X,range,order)

Y=0*X;

for i=1:length(X)
   Y(i)=legev(X(i),range,order);
end
