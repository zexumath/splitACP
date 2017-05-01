function [W,Wmat] = generatew3(opt1D,VNew,DNew,Nocc,zshift,zweight,RHSsel,sel)

Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
Ntot    = length(DNew);


if nargout ==2
    Wmat1 = zeros(NsGrid,NsGrid,select);
    Wmat2 = zeros(NsGrid,NsGrid,select);
    W     = zeros(NsGrid,select);
    for nu = 1:select
        for l = 1:Npole
            zeta1 = VNew(:,1:Nocc) * ((DNew(1:Nocc) - zshift(l)).\ (VNew(:,1:Nocc)'*RHSsel(:,nu)* hs));
            zeta2 = VNew(:,Nocc+1:Ntot) * ((DNew(Nocc+1:Ntot) - zshift(l)).\ (VNew(:,Nocc+1:Ntot)'*RHSsel(:,nu)* hs));
            for i = 1:Nocc
                fac = VNew(sel(nu),i) * zweight(l)...
                        / (zshift(l) - DNew(i));
                Wmat1(:,:,nu) = Wmat1(:,:,nu) + ...
                    (zeta1 * VNew(:,i)') * fac;
                Wmat2(:,:,nu) = Wmat2(:,:,nu) + ...
                    (zeta2 * VNew(:,i)') * fac;
            end
        end
        Wmat1(:,:,nu) = 1/(2i) * (Wmat1(:,:,nu) - Wmat1(:,:,nu)');
        Wmat2(:,:,nu) = 1/(2i) * (Wmat2(:,:,nu) - Wmat2(:,:,nu)');
        W(:,nu) = real(diag(Wmat1(:,:,nu)) + diag(Wmat2(:,:,nu)) + diag(conj(Wmat2(:,:,nu))));
    end
    Wmat = real(Wmat1 + Wmat2 + conj(Wmat2));

elseif nargout ==1
    W1  = zeros(NsGrid,select);
    W2  = zeros(NsGrid,select);
    for nu = 1:select
        for l = 1:Npole
            zeta1 = VNew(:,1:Nocc) * ((DNew(1:Nocc) - zshift(l)).\ (VNew(:,1:Nocc)'*RHSsel(:,nu)* hs));
            zeta2 = VNew(:,Nocc+1:Ntot) * ((DNew(Nocc+1:Ntot) - zshift(l)).\ (VNew(:,Nocc+1:Ntot)'*RHSsel(:,nu)* hs));
            for i = 1:Nocc
                fac = VNew(sel(nu),i) * zweight(l)...
                        / (zshift(l) - DNew(i));
                    W1(:,nu) = W1(:,nu) + ...
                        (conj(VNew(:,i)).* zeta1) * fac;
                    W2(:,nu) = W2(:,nu) + ...
                        (conj(VNew(:,i)).* zeta2) * fac;
            end
        end
    end
    W1 = 1/(2i) * (W1 - conj(W1));
    W2 = 1/(2i) * (W2 - conj(W2));
    W = W1 + W2 + conj(W2);
end