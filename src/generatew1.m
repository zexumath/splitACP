function [drhodRACPs, dPdRACPs] = generatew1(opt1D,VNew,DNew,zshift,zweight,Vcg)

Npole   = opt1D.Npole;
hs      = opt1D.hs;
NsGrid  = opt1D.NsGrid;
NsCell  = opt1D.atom.NsCell;
dPdRACPs    = zeros(NsGrid,NsGrid,size(Vcg,2));
Ntot    = length(DNew);

if nargout ==2
    dPdRACPs    = zeros(NsGrid,NsGrid,size(Vcg,2));
    drhodRACPs  = zeros(NsGrid, size(Vcg,2));
    for j = 1:NsCell
        res3 = zeros(NsGrid, NsGrid);
        for l = 1:Npole
            Gl = @(x) VNew *(diag(DNew - zshift(l))\(VNew'*x))* hs;
            tmp3 = zeros(NsGrid, NsGrid);
            for i = 1:Ntot
                tmp2 = Vcg(:,j).*VNew(:,i);
                tmp3 = tmp3 + Gl( tmp2 ) / (zshift(l) - DNew(i)) *VNew(:,i)';
            end
            res3 = res3 + zweight(l) * tmp3;
        end
        res3 = 1/(2i)*(res3-res3');
        dPdRACPs(:,:,j) = res3;
        drhodRACPs (:,j) = diag(drhodRACPs(:,:,j));
    end

elseif nargout ==1
    drhodRACPs  = zeros(NsGrid, size(Vcg,2));
    for j = 1:NsCell
        res1 = zeros(NsGrid, 1);
        for l = 1:Npole
            Gl = @(x) VNew *(diag(DNew - zshift(l))\(VNew'*x))* hs;
            tmp1 = zeros(NsGrid, 1);
            for i = 1:Ntot
                tmp2 = Vcg(:,j).*VNew(:,i);
                tmp1 = tmp1 + Gl( tmp2 ) / (zshift(l) - DNew(i)) .*conj(VNew(:,i));
            end
            res1 = res1 + zweight(l) * tmp1;
        end
        res1 = 1/(2i)*(res1-conj(res1));
        drhodRACPs (:,j) = res1;
    end
end