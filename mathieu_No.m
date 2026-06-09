function [R, dRdu, coeffs] = mathieu_No(n, u, q, N)
%MATHIEU_NO  Radial odd Mathieu function (2nd kind) No_n(u,q).
%   R=MATHIEU_NO(n,u,q): Bessel Y product series, same as Jo but J->Y.
%   No_n is the outgoing-wave complement to Jo_n. Diverges at u→0.
validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
u=u(:)'; q=abs(q);
if n==0, R=zeros(size(u)); dRdu=zeros(size(u)); coeffs=zeros(0,1); return; end
if any(u<=0), warning('MATHIEU_NO:SmallU','No_n diverges for u→0.'); end
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
        Yr1=bessely(r,v1); Yr1_1=bessely(r+1,v1);
        Yr2=bessely(r,v2); Yr2_1=bessely(r+1,v2);
        R=R+sf*Br.*(Yr1.*Yr2_1-Yr1_1.*Yr2);
        if nargout>=2
            Yr1p=bessel_deriv('Y',r,v1); Yr1_1p=bessel_deriv('Y',r+1,v1);
            Yr2p=bessel_deriv('Y',r,v2); Yr2_1p=bessel_deriv('Y',r+1,v2);
            dRdu=dRdu+sf*Br.*(-v1.*Yr1p.*Yr2_1+v2.*Yr1.*Yr2_1p+v1.*Yr1_1p.*Yr2-v2.*Yr1_1.*Yr2p);
        end
    end
    B1=coeffs(1); if abs(B1)>1e-15, R=R/B1; if nargout>=2, dRdu=dRdu/B1; end, end
else
    m=(n-2)/2;
    for r=0:(M-1)
        Br=coeffs(r+1); if abs(Br)<1e-15, continue; end
        sf=(-1)^(r+m);
        Yr1=bessely(r,v1); Yr1_2=bessely(r+2,v1);
        Yr2=bessely(r,v2); Yr2_2=bessely(r+2,v2);
        R=R+sf*Br.*(Yr1.*Yr2_2-Yr1_2.*Yr2);
        if nargout>=2
            Yr1p=bessel_deriv('Y',r,v1); Yr1_2p=bessel_deriv('Y',r+2,v1);
            Yr2p=bessel_deriv('Y',r,v2); Yr2_2p=bessel_deriv('Y',r+2,v2);
            dRdu=dRdu+sf*Br.*(-v1.*Yr1p.*Yr2_2+v2.*Yr1.*Yr2_2p+v1.*Yr1_2p.*Yr2-v2.*Yr1_2.*Yr2p);
        end
    end
    B2=coeffs(1); if abs(B2)>1e-15, R=R/B2; if nargout>=2, dRdu=dRdu/B2; end, end
end
end
