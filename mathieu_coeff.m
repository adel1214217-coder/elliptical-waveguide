function [coeffs, char_val] = mathieu_coeff(n, q, func_type, N)
%MATHIEU_COEFF  Fourier expansion coefficients for angular Mathieu functions.
%   coeffs = MATHIEU_COEFF(n,q,'ce') returns A_r for ce_n(v,q)=Σ A_k cos(k·v).
%   coeffs = MATHIEU_COEFF(n,q,'se') returns B_r for se_n(v,q)=Σ B_k sin(k·v).
%   [coeffs,a] = MATHIEU_COEFF(...) also returns characteristic value a_n or b_n.
%
%   Normalization (NIST DLMF 28.2): ∫ce_n²=π (n>0), 2π (n=0); ∫se_n²=π.
%   Sign: ce_{2m}(0)>0, ce_{2m+1}'(0)>0, se_n'(0)>0.
%
%   Method: Solve tridiagonal eigenvalue problem from Fourier expansion
%   of angular Mathieu equation d²y/dv²+(a-2q·cos2v)y=0.
%   Input: n=order, q=parameter, func_type='ce'|'se', N=truncation (auto).
%   Output: coeffs (N×1), char_val (a_n or b_n).

validateattributes(n,{'numeric'},{'integer','nonnegative','scalar'});
validateattributes(q,{'numeric'},{'real','nonnegative','scalar'});
func_type = validatestring(func_type, {'ce','se'});

if strcmp(func_type,'se') && n==0
    coeffs=zeros(0,1); char_val=NaN; return;
end

if nargin<4, N=max([40, n+30, ceil(6*sqrt(q)+25)]); end
q=abs(q);

% Build matrix, solve eigenproblem
if strcmp(func_type,'ce'), T=tridiag_ce(n,q,N); else, T=tridiag_se(n,q,N); end
[V,D]=eig(T); [evals,idx]=sort(diag(D));

m=floor(n/2);
if m+1>length(evals)
    error('MATHIEU_COEFF:TruncationSmall',...
          'N=%d too small for n=%d, q=%.4f. Increase N.',N,n,q);
end
char_val=evals(m+1); coeffs_raw=V(:,idx(m+1));

% Normalization
if strcmp(func_type,'ce')
    if n==0, int_val=2*pi*coeffs_raw(1)^2+pi*sum(coeffs_raw(2:end).^2); target=2*pi;
    elseif mod(n,2)==0, int_val=2*pi*coeffs_raw(1)^2+pi*sum(coeffs_raw(2:end).^2); target=pi;
    else, int_val=pi*sum(coeffs_raw.^2); target=pi;
    end
    % Sign: ce_{2m}(0)>0 => sum(coeffs)>0 for even n
    if mod(n,2)==0
        if sum(coeffs_raw)<0, coeffs_raw=-coeffs_raw; end
    else
        r=(0:N-1)'; kvals=2*r+1;
        if -sum(kvals.*coeffs_raw.*sin(kvals*1e-6))<0, coeffs_raw=-coeffs_raw; end
    end
else
    int_val=pi*sum(coeffs_raw.^2); target=pi;
    % Sign: se_n'(0)>0
    r=(0:N-1)';
    if mod(n,2)==1, deriv=sum((2*r+1).*coeffs_raw);
    else, deriv=sum((2*r+2).*coeffs_raw);
    end
    if deriv<0, coeffs_raw=-coeffs_raw; end
end
coeffs=coeffs_raw*sqrt(target/int_val);
end
