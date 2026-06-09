function [y, coeffs, a] = mathieu_ce(n, v, q, N)
%MATHIEU_CE  Angular even Mathieu function ce_n(v,q).
%   y=MATHIEU_CE(n,v,q): n=0,1,2,..., v=radians (scalar/vector), q>=0.
%   [y,coeffs,a]=MATHIEU_CE(...) also returns Fourier coefficients and a_n(q).
%   n even: ce_n=Σ A_{2r} cos(2r·v); n odd: ce_n=Σ A_{2r+1} cos((2r+1)·v).
%   Vectorized: y=coeffs'*cos(k_vals*v). Normalized per NIST DLMF 28.2.
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
v=v(:)';
if nargin<4, N=max([40,n+30,ceil(6*sqrt(q)+25)]); end
[coeffs,a]=mathieu_coeff(n,q,'ce',N); M=length(coeffs);
if mod(n,2)==0, kv=(0:2:(2*M-2))'; else, kv=(1:2:(2*M-1))'; end
y=coeffs.'*cos(kv*v);
end
