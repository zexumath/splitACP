function [drhodRACPb, drhodRACPb1, drhodRACPb2] = generatew2(opt1D,VNew,DNew,Nocc,zshift,zweight,gpfunc1,RHSsel,sel)

Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
NsCell  = opt1D.atom.NsCell;
Ntot    = length(DNew);

if(0)
    dPdRACPb1  = zeros(NsGrid,NsGrid,NsCell);
    dPdRACPb2  = zeros(NsGrid,NsGrid,NsCell);
    for j = 1:NsCell
        for nu = 1:select
            for l = 1:Npole
                zeta1 = VNew(:,1:Nocc) * ((DNew(1:Nocc) - zshift(l)).\ (VNew(:,1:Nocc)'*RHSsel(:,nu)* hs));
                zeta2 = VNew(:,Nocc+1:Ntot) * ((DNew(Nocc+1:Ntot) - zshift(l)).\ (VNew(:,Nocc+1:Ntot)'*RHSsel(:,nu)* hs));
                for i = 1:Nocc
                    fac = gpfunc1(VNew(:,i),j,sel(nu)) * zweight(l)...
                        / (zshift(l) - DNew(i));
                    dPdRACPb1(:,j) = dPdRACPb1(:,j) + ...
                        ( zeta1*VNew(:,i)') * fac;
                    dPdRACPb2(:,j) = dPdRACPb2(:,j) + ...
                        ( zeta2*VNew(:,i)') * fac;
                end
            end
        end
        dPdRACPb1(:,j) = 1/(2i)* (dPdRACPb1(:,j) -(dPdRACPb1(:,j))');
        dPdRACPb2(:,j) = 1/(2i)* (dPdRACPb2(:,j) -(dPdRACPb2(:,j))');
    end
    dPdRACPb = dPdRACPb1 + dPdRACPb2 + (dPdRACPb2)';
    
else
    drhodRACPb1  = zeros(NsGrid,NsCell);
    drhodRACPb2  = zeros(NsGrid,NsCell);
    for j = 1:NsCell
        for nu = 1:select
            for l = 1:Npole
                zeta1 = VNew(:,1:Nocc) * ((DNew(1:Nocc) - zshift(l)).\ (VNew(:,1:Nocc)'*RHSsel(:,nu)* hs));
                zeta2 = VNew(:,Nocc+1:Ntot) * ((DNew(Nocc+1:Ntot) - zshift(l)).\ (VNew(:,Nocc+1:Ntot)'*RHSsel(:,nu)* hs));
                for i = 1:Nocc
                    fac = gpfunc1(VNew(:,i),j,sel(nu)) * zweight(l)...
                        / (zshift(l) - DNew(i));
                    drhodRACPb1(:,j) = drhodRACPb1(:,j) + ...
                        (conj(VNew(:,i)).* zeta1) * fac;
                    drhodRACPb2(:,j) = drhodRACPb2(:,j) + ...
                        (conj(VNew(:,i)).* zeta2) * fac;
                end
            end
        end
        drhodRACPb1(:,j) = 1/(2i)* (drhodRACPb1(:,j) -conj(drhodRACPb1(:,j)));
        drhodRACPb2(:,j) = 1/(2i)* (drhodRACPb2(:,j) -conj(drhodRACPb2(:,j)));
    end
    drhodRACPb = drhodRACPb1 + drhodRACPb2 + conj(drhodRACPb2);
end

end
