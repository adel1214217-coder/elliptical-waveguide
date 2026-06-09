function [R, dRdu, coeffs] = mathieu_Jo(n, u, q, N)
%MATHIEU_JO  Radial odd Mathieu function (1st kind) Jo_n(u,q).
%   R=MATHIEU_JO(n,u,q): Bessel J product series (NIST DLMF 28.20).
%   n odd: Jo=(1/B1)Σ(-1)^{r+m}B_{2r+1}[J_r(v1)J_{r+1}(v2)-J_{r+1}(v1)J_r(v2)].
%   n even: Jo=(1/B2)Σ(-1)^{r+m}B_{2r+2}[J_r(v1)J_{r+2}(v2)-J_{r+2}(v1)J_r(v2)].
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
u=u(:)'; q=abs(q);
if n==0, R=zeros(size(u)); dRdu=zeros(size(u)); coeffs=zeros(0,1); return; end
ang_N=max([50,n+40,ceil(6*sqrt(q)+35)]);
[coeffs,~]=mathieu_coeff(n,q,'se',ang_N); M=length(coeffs);
v1=sqrt(q)*exp(-u); v2=sqrt(q)*exp(u);
if nargin<4, N=max(40,ceil(sqrt(q)*exp(max(u))+20)); end
R=zeros(size(u)); dRdu=zeros(size(u));
if mod(n,2)==1
    m=(n-1)/2;
    for r=0:(M-1)
        Br=coeffs(r+1); if abs(Br)<1e-15, continue; end
        sf=(-1)^(r+m);
        Jr1=besselj(r,v1); Jr1_1=besselj(r+1,v1);
        Jr2=besselj(r,v2); Jr2_1=besselj(r+1,v2);
        R=R+sf*Br.*(Jr1.*Jr2_1-Jr1_1.*Jr2);
        if nargout>=2
            Jr1p=bessel_deriv('J',r,v1); Jr1_1p=bessel_deriv('J',r+1,v1);
            Jr2p=bessel_deriv('J',r,v2); Jr2_1p=bessel_deriv('J',r+1,v2);
            dRdu=dRdu+sf*Br.*(-v1.*Jr1p.*Jr2_1+v2.*Jr1.*Jr2_1p+v1.*Jr1_1p.*Jr2-v2.*Jr1_1.*Jr2p);
        end
    end
    B1=coeffs(1); if abs(B1)>1e-15, R=R/B1; if nargout>=2, dRdu=dRdu/B1; end, end
else
    m=(n-2)/2;
    for r=0:(M-1)
        Br=coeffs(r+1); if abs(Br)<1e-15, continue; end
        sf=(-1)^(r+m);
        Jr1=besselj(r,v1); Jr1_2=besselj(r+2,v1);
        Jr2=besselj(r,v2); Jr2_2=besselj(r+2,v2);
        R=R+sf*Br.*(Jr1.*Jr2_2-Jr1_2.*Jr2);
        if nargout>=2
            Jr1p=bessel_deriv('J',r,v1); Jr1_2p=bessel_deriv('J',r+2,v1);
            Jr2p=bessel_deriv('J',r,v2); Jr2_2p=bessel_deriv('J',r+2,v2);
            dRdu=dRdu+sf*Br.*(-v1.*Jr1p.*Jr2_2+v2.*Jr1.*Jr2_2p+v1.*Jr1_2p.*Jr2-v2.*Jr1_2.*Jr2p);
        end
    end
    B2=coeffs(1); if abs(B2)>1e-15, R=R/B2; if nargout>=2, dRdu=dRdu/B2; end, end
end
end
