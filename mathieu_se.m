function [y, coeffs, b] = mathieu_se(n, v, q, N)
%MATHIEU_SE  Angular odd Mathieu function se_n(v,q).
%   y=MATHIEU_SE(n,v,q): n=1,2,3,..., v=radians, q>=0. (n=0→y=0).
%   [y,coeffs,b]=MATHIEU_SE(...) also returns coefficients and b_n(q).
%   n odd: se_n=Σ B_{2r+1} sin((2r+1)·v); n even: se_n=Σ B_{2r+2} sin((2r+2)·v).
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
v=v(:)';
if n==0, y=zeros(size(v)); coeffs=zeros(0,1); b=NaN; return; end
if nargin<4, N=max([40,n+30,ceil(6*sqrt(q)+25)]); end
[coeffs,b]=mathieu_coeff(n,q,'se',N); M=length(coeffs);
if mod(n,2)==1, kv=(1:2:(2*M-1))'; else, kv=(2:2:(2*M))'; end
y=coeffs.'*sin(kv*v);
end
