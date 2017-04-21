function chi0s = applychi0s2(opt1D,umat, DNew,VNew,efermi,Npole)

T       = opt1D.temperature;
Gap     = 0.0;
DeltaE  = max(2.0, DNew(end) - DNew(1));
mu      = efermi;
hs      = opt1D.hs;
NsGrid  = opt1D.NsGrid;


[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);

Ntot    = length(DNew);

[~,m] = size(umat);

chi0s = zeros(NsGrid,m);

for ncol = 1:m
    res = zeros(NsGrid, 1);
    for l = 1:Npole
        Gl = VNew / diag(DNew - zshift(l)) *VNew' * hs;
        tmp1 = zeros(NsGrid, 1);
        for i = 1:Ntot
            tmp2 = umat(:,ncol).* VNew(:,i);
            tmp1 = tmp1 + Gl * tmp2 / (zshift(l) - DNew(i)) .* conj(VNew(:,i));
        end
        res = res + zweight(l) * tmp1;
    end
    chi0s(:,ncol) = 1/(2i)*( res - conj(res) );
end


end