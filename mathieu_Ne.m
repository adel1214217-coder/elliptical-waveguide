function [R, dRdu, coeffs] = mathieu_Ne(n, u, q, N)
%MATHIEU_NE  Radial even Mathieu function (2nd kind) Ne_n(u,q).
%   R=MATHIEU_NE(n,u,q): Bessel Y product series, same as Je but J->Y.
%   Ne_n is the outgoing-wave complement to Je_n. Diverges at u→0.
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
u=u(:)'; q=abs(q);
if any(u<=0), warning('MATHIEU_NE:SmallU','Ne_n diverges for u→0.'); end
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
        Yr1=bessely(r,v1); Yr2=bessely(r,v2);
        R=R+sf*Ar*Yr1.*Yr2;
        if nargout>=2
            Yr1p=bessel_deriv('Y',r,v1); Yr2p=bessel_deriv('Y',r,v2);
            dRdu=dRdu+sf*Ar.*(-v1.*Yr1p.*Yr2+v2.*Yr1.*Yr2p);
        end
    end
    A0=coeffs(1); if abs(A0)>1e-15, R=R/A0; if nargout>=2, dRdu=dRdu/A0; end, end
else
    m=(n-1)/2;
    for r=0:(M-1)
        Ar=coeffs(r+1); if abs(Ar)<1e-15, continue; end
        sf=(-1)^(r+m);
        Yr1=bessely(r,v1); Yr1_1=bessely(r+1,v1);
        Yr2=bessely(r,v2); Yr2_1=bessely(r+1,v2);
        R=R+sf*Ar.*(Yr1.*Yr2_1+Yr1_1.*Yr2);
        if nargout>=2
            Yr1p=bessel_deriv('Y',r,v1); Yr1_1p=bessel_deriv('Y',r+1,v1);
            Yr2p=bessel_deriv('Y',r,v2); Yr2_1p=bessel_deriv('Y',r+1,v2);
            dRdu=dRdu+sf*Ar.*(-v1.*Yr1p.*Yr2_1+v2.*Yr1.*Yr2_1p-v1.*Yr1_1p.*Yr2+v2.*Yr1_1.*Yr2p);
        end
    end
    A1=coeffs(1); if abs(A1)>1e-15, R=R/A1; if nargout>=2, dRdu=dRdu/A1; end, end
end
end
