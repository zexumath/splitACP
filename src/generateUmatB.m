function [dPdRbs, drhodRbs] = generateUmatB(opt1D,VNew,DNew,gpfunc,zshift,zweight,RHSsel,sel)

Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
Ntot    = length(DNew);
NsCell  = opt1D.atom.NsCell;
dPdRbs  = zeros(NsGrid,NsGrid,NsCell);
drhodRbs= zeros(NsGrid,NsCell);

G = zeros(NsGrid,NsCell,Ntot);
for na = 1:NsCell
    for i = 1:Ntot
        G(:,na,i) = gpfunc(VNew(:,i),na);
    end
end

for na = 1:NsCell
    for nu = 1:select
        for l = 1:Npole
            zeta = VNew / diag(DNew - zshift(l)) * (VNew'*RHSsel(:,nu)* hs);
            for i = 1:Ntot
                dPdRbs(:,:,na) = dPdRbs(:,:,na) + zweight(l)/ (zshift(l) - DNew(i)) ...
                    * zeta * VNew(:,i)'* G(sel(nu),na,i);
            end
        end
    end
    dPdRbs(:,:,na) = 1/(2i) * (dPdRbs(:,:,na) - dPdRbs(:,:,na)');
    fprintf('na = %d\n',na);
    drhodRbs(:,na) = real(diag(dPdRbs(:,:,na)));
end
end