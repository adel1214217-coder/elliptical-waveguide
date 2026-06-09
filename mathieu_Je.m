function [R, dRdu, coeffs] = mathieu_Je(n, u, q, N)
%MATHIEU_JE  Radial even Mathieu function (1st kind) Je_n(u,q).
%   R=MATHIEU_JE(n,u,q): Bessel J product series (NIST DLMF 28.20).
%   v1=√q·e^{-u}, v2=√q·e^{u}. n even: Je=(1/A0)Σ(-1)^{r+m}A_{2r}J_r(v1)J_r(v2).
%   n odd: Je=(1/A1)Σ(-1)^{r+m}A_{2r+1}[J_r(v1)J_{r+1}(v2)+J_{r+1}(v1)J_r(v2)].
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
u=u(:)'; q=abs(q);
ang_N=max([50,n+40,ceil(6*sqrt(q)+35)]);
[coeffs,~]=mathieu_coeff(n,q,'ce',ang_N); M=length(coeffs);
v1=sqrt(q)*exp(-u); v2=sqrt(q)*exp(u);
if nargin<4, N=max(40,ceil(sqrt(q)*exp(max(u))+20)); end
R=zeros(size(u)); dRdu=zeros(size(u));
if mod(n,2)==0
    m=n/2;
    for r=0:(M-1)
        Ar=coeffs(r+1); if abs(Ar)<1e-15, continue; end
        sf=(-1)^(r+m);
        Jr1=besselj(r,v1); Jr2=besselj(r,v2);
        R=R+sf*Ar*Jr1.*Jr2;
        if nargout>=2
            Jr1p=bessel_deriv('J',r,v1); Jr2p=bessel_deriv('J',r,v2);
            dRdu=dRdu+sf*Ar.*(-v1.*Jr1p.*Jr2+v2.*Jr1.*Jr2p);
        end
    end
    A0=coeffs(1); if abs(A0)>1e-15, R=R/A0; if nargout>=2, dRdu=dRdu/A0; end, end
else
    m=(n-1)/2;
    for r=0:(M-1)
        Ar=coeffs(r+1); if abs(Ar)<1e-15, continue; end
        sf=(-1)^(r+m);
        Jr1=besselj(r,v1); Jr1_1=besselj(r+1,v1);
        Jr2=besselj(r,v2); Jr2_1=besselj(r+1,v2);
        R=R+sf*Ar.*(Jr1.*Jr2_1+Jr1_1.*Jr2);
        if nargout>=2
            Jr1p=bessel_deriv('J',r,v1); Jr1_1p=bessel_deriv('J',r+1,v1);
            Jr2p=bessel_deriv('J',r,v2); Jr2_1p=bessel_deriv('J',r+1,v2);
            dRdu=dRdu+sf*Ar.*(-v1.*Jr1p.*Jr2_1+v2.*Jr1.*Jr2_1p-v1.*Jr1_1p.*Jr2+v2.*Jr1_1.*Jr2p);
        end
    end
    A1=coeffs(1); if abs(A1)>1e-15, R=R/A1; if nargout>=2, dRdu=dRdu/A1; end, end
end
end
