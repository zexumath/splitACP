function [drhodRACPb] = generatew2(opt1D,bV,bxV,VNew,DNew,Nocc,zshift,zweight,zetatilde,sel)

Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
NsCell  = opt1D.atom.NsCell;
Ntot    = length(DNew);

b       = opt1D.nlpp;
bx      = opt1D.nlppx;

% if(nargout==2)
%     dPdRACPb1  = zeros(NsGrid,NsGrid,NsCell);
%     dPdRACPb2  = zeros(NsGrid,NsGrid,NsCell);
%     for j = 1:NsCell
%         for nu = 1:select
%             for l = 1:Npole
%                 zeta1 = VNew(:,1:Nocc) * zetatilde(1:Nocc,l,nu);
%                 zeta2 = VNew(:,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,l,nu);
%                 for i = 1:Nocc
%                     fac = (VNew(sel(nu),i).*gp(sel(nu),j) ...
%                         +(-1)*b(sel(nu),j)*bxV(j,i) + (-1)*bx(sel(nu),j)*bV(j,i))...
%                         * zweight(l)/ (zshift(l) - DNew(i));
%                     dPdRACPb1(:,j) = dPdRACPb1(:,j) + ...
%                         ( zeta1*VNew(:,i)') * fac;
%                     dPdRACPb2(:,j) = dPdRACPb2(:,j) + ...
%                         ( zeta2*VNew(:,i)') * fac;
%                 end
%             end
%         end
%         dPdRACPb1(:,j) = 1/(2i)* (dPdRACPb1(:,j) -(dPdRACPb1(:,j))');
%         dPdRACPb2(:,j) = 1/(2i)* (dPdRACPb2(:,j) -(dPdRACPb2(:,j))');
%     end
%     dPdRACPb = dPdRACPb1 + dPdRACPb2 + (dPdRACPb2)';
%
% else
%     drhodRACPb1  = zeros(NsGrid,NsCell);
%     drhodRACPb2  = zeros(NsGrid,NsCell);
%     for j = 1:NsCell
%         for nu = 1:select
%             for l = 1:Npole
%                 zeta1 = VNew(:,1:Nocc) * zetatilde(1:Nocc,l,nu);
%                 zeta2 = VNew(:,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,l,nu);
%                 for i = 1:Nocc
%                     fac = (VNew(sel(nu),i).*gp(sel(nu),j) ...
%                         +(-1)*b(sel(nu),j)*bxV(j,i) + (-1)*bx(sel(nu),j)*bV(j,i))...
%                         * zweight(l)/ (zshift(l) - DNew(i));
%                     drhodRACPb1(:,j) = drhodRACPb1(:,j) + ...
%                         (conj(VNew(:,i)).* zeta1) * fac;
%                     drhodRACPb2(:,j) = drhodRACPb2(:,j) + ...
%                         (conj(VNew(:,i)).* zeta2) * fac;
%                 end
%             end
%         end
%         drhodRACPb1(:,j) = 1/(2i)* (drhodRACPb1(:,j) -conj(drhodRACPb1(:,j)));
%         drhodRACPb2(:,j) = 1/(2i)* (drhodRACPb2(:,j) -conj(drhodRACPb2(:,j)));
%     end
%     drhodRACPb = drhodRACPb1 + drhodRACPb2 + conj(drhodRACPb2);
% end
drhodRACPb1  = zeros(NsGrid,NsCell);
drhodRACPb2  = zeros(NsGrid,NsCell);
for j = 1:NsCell
    for nu = 1:select
        if b(sel(nu),j)==0 && b(sel(nu),j)==0
            continue
        else
            for l = 1:Npole
                zeta1 = VNew(:,1:Nocc) * zetatilde(1:Nocc,nu,l);
                zeta2 = VNew(:,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,nu,l);
                for i = 1:Nocc
                    fac = ((-1)*b(sel(nu),j)*bxV(j,i) + (-1)*bx(sel(nu),j)*bV(j,i))...
                        * zweight(l)/ (zshift(l) - DNew(i));
                    drhodRACPb1(:,j) = drhodRACPb1(:,j) + ...
                        (conj(VNew(:,i)).* zeta1) * fac;
                    drhodRACPb2(:,j) = drhodRACPb2(:,j) + ...
                        (conj(VNew(:,i)).* zeta2) * fac;
                end
            end
        end
    end
    drhodRACPb1(:,j) = 1/(2i)* (drhodRACPb1(:,j) -conj(drhodRACPb1(:,j)));
    drhodRACPb2(:,j) = 1/(2i)* (drhodRACPb2(:,j) -conj(drhodRACPb2(:,j)));
end
drhodRACPb = drhodRACPb1 + drhodRACPb2 + conj(drhodRACPb2);

end
