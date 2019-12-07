function X=ProjectOnSn(X)
% Project a point or a set of points onto Sn.
%
% INPUT:
%   - X : N-by-d array of point coordinates
%   
% OUTPUT:
%   - X : N-by-d array of point coordinates such that norm(X(i,:))=1
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


X=bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));


