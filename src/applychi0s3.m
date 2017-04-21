function chi0smat = applychi0s3(opt1D,umat, DNew,VNew,efermi,Npole)

T       = opt1D.temperature;
Gap     = 0.0;
DeltaE  = max(2.0, DNew(end) - DNew(1));
mu      = efermi;
hs      = opt1D.hs;
NsGrid  = opt1D.NsGrid;


[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);

Nt      = length(DNew);

[~,m] = size(umat);

chi0smat = zeros(NsGrid,NsGrid,m);

for ncol = 1:m
    res = zeros(NsGrid, NsGrid);
    for l = 1:Npole
        Gl = @(x) VNew *( diag(DNew - zshift(l))\ (VNew'*x)) * hs;
        tmp1 = zeros(NsGrid, NsGrid);
        for i = 1:Nt
            tmp2 = umat(:,ncol).* VNew(:,i);
            tmp1 = tmp1 + Gl( tmp2) / (zshift(l) - DNew(i)) * VNew(:,i)';
        end
        res = res + zweight(l) * tmp1;
    end
    chi0smat(:,:,ncol) = 1/(2i)*( res - res' );
end


end