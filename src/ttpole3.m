% This script test on the integral of f(z) (z-H1)^-1 g (z-H2)^-1
% See if the actual results matches with integral, either with or without
% the compression.

Npole   = opt1D.Npole;
T       = opt1D.temperature;
Gap     = 0.0;
DeltaE  = max(2.0, DNew(end) - DNew(1));
mu      = efermi;

% direct sum
result1 = zeros(NsGrid,NsCell);
for J =1:NsCell
    for k = 1:Nocc
        for j = Nocc+1:Ntot
            result1(:,J) = result1(:,J) + (occ(j) - occ(k)) / (DNew(j) - DNew(k)) ...
                * VNew(:,j).*conj(VNew(:,k)) * (hs*VNew(:,j)'*gpfunc(VNew(:,k),J));
        end
    end
end

% integral with pole expansion: error 1e-10
[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);
result2 = zeros(NsGrid,NsCell);
for j = 1:NsCell
    res = zeros(NsGrid, 1);
    for l = 1:Npole
        Gl = @(x) VNew(:,Nocc+1:Ntot) *((DNew(Nocc+1:Ntot) - zshift(l)).\(VNew(:,Nocc+1:Ntot)'*x))* hs;
        tmp1 = zeros(NsGrid, 1);
        for i = 1:Nocc
            tmp2 = gpfunc(VNew(:,i),j);
            tmp1 = tmp1 + Gl( tmp2 ) / (zshift(l) - DNew(i)) .* conj(VNew(:,i));
        end
        res = res + zweight(l) * tmp1;
    end
    res = 1/(2i)*(res-conj(res));
    result2(:,j) = res;
end




