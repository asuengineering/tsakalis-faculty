%function Y=legev_bvec(X,range,order_index)
% X: Vector of elements
% range: Vector 2 elements
% order_index: vector with degrees of polynomial, 
% Y: matrix with k-columns (k=length(order_index), the evaluation of Chevyshev polynomial 
%    of degre order_index(i)
%
%i.e., P(x)=C0(x)+C1(x)+C2(x)  
%	order_index=[0;1;2]
% where C0(x) is Chevyshev polynomial degree 0
%		 C1(x) is Chevyshev polynomial degree 1
%       C2(x) is Chevyshev polynomial degree 2

function Y=legev_bv(X,range,ord_indx)


order=length(ord_indx);


Y=zeros(length(X),order);

for i=1:length(X)
   for k=1:order
      Y(i,k)=legev(X(i),range,ord_indx(k));
   end
end
