function res3 = getdbPb(opt1D,bV,bxV,I,J,DNew,VNew,zshift,zweight,Vcg,gpfunc)
Npole   = opt1D.Npole;
hs      = opt1D.hs;
NsGrid  = opt1D.NsGrid;
NsCell  = opt1D.atom.NsCell;
Ntot    = length(DNew);
res3    = 0;
for j = 1:NsCell
    
    for l = 1:Npole
        bGl  = @(x) bV(I,:) *((DNew - zshift(l)).\(VNew'*x))* hs;
        bxGl = @(x) bxV(I,:) *((DNew - zshift(l)).\(VNew'*x))* hs;
        tmp3 = 0;
        for i = 1:Ntot
            tmp2 = Vcg(:,J).*VNew(:,i) + gpfunc(VNew(:,i),J);
            tmp3 = tmp3 + bxGl( tmp2 ) / (zshift(l) - DNew(i)) *bV(I,i)';
            tmp3 = tmp3 + bGl( tmp2 )  / (zshift(l) - DNew(i)) *bxV(I,i)';
        end
        res3 = res3 + zweight(l) * tmp3;
    end
    res3 = 1/(2i)*(res3-res3');
end

end
