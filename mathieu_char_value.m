function [a, b] = mathieu_char_value(n, q, N)
%MATHIEU_CHAR_VALUE  Characteristic values a_n(q), b_n(q) of Mathieu equation.
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
if nargin<3, N=max([40,n+30,ceil(6*sqrt(q)+25)]); end
T=tridiag_ce(n,q,N); evals=sort(eig(T)); m=floor(n/2)+1; a=evals(min(m,length(evals)));
if n==0, b=NaN; else, T=tridiag_se(n,q,N); evals=sort(eig(T)); b=evals(min(m,length(evals))); end
end
