function chi0s = applychi0s1(opt1D,g, DNew,VNew,occ,efermi)

hs      = opt1D.hs;
NsGrid  = opt1D.NsGrid;
Tbeta   = opt1D.temperature*3.166815d-6;
Ntot    = length(DNew);
Ne      = opt1D.Ne;
gpfunc  = @(x,i) g(:,i).*x;

ita     = zeros(NsGrid,1);
theta   = zeros(NsGrid,1);
na      = Ne;

chi0s = 0;
for i = 1:Ntot
    for j = 1:Ntot
        if(abs(occ(i) - occ(j))>1e-8 )
            fac = (occ(i) - occ(j))/(DNew(i) - DNew(j));
        else
            fac = fermidiracderiv(DNew(i),efermi,Tbeta);
        end
        chi0s = chi0s + fac * ( hs * VNew(:,j)'*(g.*VNew(:,i)) )...
            *VNew(:,j) .*conj(VNew(:,i));
    end
end

end