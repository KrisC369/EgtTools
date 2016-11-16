function Out = solveN(A)
syms x
[~,agents]=size(A);
nchoosekCache = nan(agents+1,agents+1);
n=agents/2;
F = zeros(n+1,2);
L=A;
F(1,1) = L(1);
F(n+1,2) = L(end);
for i = 2:((2*n)-1)
   F(fix(i/2)+1,mod(i,2)+1) = L(i); 
end

ue1 = 0;

for i=0:n-1
    if(~isnan(nchoosekCache(n-1,i+1)))
        fac = nchoosekCache(n-1,i+1);
    else
        fac = nchoosek(n-1,i);
        nchoosekCache(n-1,i+1) = fac;
    end
    fval = F(i+1,1);
    ue1 = ue1+ ((x^(n-i-1))*((1-x)^(i))*fac*fval);
end

uexpart1 = 0;
uexpart2 = 0;

for i=1:2
   uexpart1 = uexpart1 + (x^(n-((i-1) * n))) * (1-x)^((i-1) * n) * F(1+(i-1)*n,i);
end

for j = 1:n-1
    coef = 0;
    for k=0:1
        if(~isnan(nchoosekCache(n-1,j-k+1)))
            fac = nchoosekCache(n-1,j-k+1);
        else
            fac = nchoosek(n-1,j-k);
            nchoosekCache(n-1,j-k+1) = fac;
        end
       coef = coef + F(j+1,k+1) * fac;
    end
    uexpart2 = uexpart2+((x^(n-j))*((1-x)^j)*coef);
end
uexx = uexpart1+uexpart2;
eqn = x*(ue1 - uexx);
T = vpasolve(eqn, x);

tol=10^-5;
T=T(T==real(T));
T = T(T<1&T>0);
if(~isempty(T))
    T = T(abs(T-1)>=tol);
end
Out = double(T);
end