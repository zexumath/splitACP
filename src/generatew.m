function [Ws, Wmat] = generatew(opt1D,VNew,DNew,zshift,zweight,RHSsel,sel)

Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
Ws      = zeros(NsGrid,select);
Ntot    = length(DNew);

if nargout ==2
    Wmat = zeros(NsGrid,NsGrid,select);
    for nu = 1:select
        for l = 1:Npole
            zeta = VNew *( (DNew - zshift(l)).\ (VNew'*RHSsel(:,nu)* hs) );
            for i = 1:Ntot
                Wmat(:,:,nu) = Wmat(:,:,nu) + ...
                    zeta * VNew(:,i)' * VNew(sel(nu),i) * zweight(l)...
                    / (zshift(l) - DNew(i));
            end
        end
        Wmat(:,:,nu) = 1/(2i) * (Wmat(:,:,nu) - Wmat(:,:,nu)');
        Ws(:,nu) = real(diag(Wmat(:,:,nu)));
    end
elseif nargout ==1
    Ws = zeros(NsGrid,select);
    for nu = 1:select
        for l = 1:Npole
            zeta = VNew *( (DNew - zshift(l)).\ (VNew'*RHSsel(:,nu)* hs) );
            for i = 1:Ntot
                Ws(:,nu) = Ws(:,nu) + ...
                    (conj(VNew(:,i)).* zeta) * VNew(sel(nu),i) * zweight(l)...
                    / (zshift(l) - DNew(i));
            end
        end
    end
    Ws = 1/(2i)* (Ws -conj(Ws));
end
end