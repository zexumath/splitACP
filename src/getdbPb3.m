function W = getdbPb3(opt1D,bV,bxV,DNew,VNew,Nocc,zshift,zweight,bzetatilde, bxzetatilde,sel)
Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
Ntot    = length(DNew);
NsCell  = opt1D.atom.NsCell;

W1 = zeros(NsCell,select);
W2 = zeros(NsCell,select);
W  = zeros(NsCell,select);

for nu = 1:select
    for l = 1:Npole
        bzeta1  = bzetatilde( 1,:,l,nu);
        bxzeta1 = bxzetatilde(1,:,l,nu);
        bzeta2  = bzetatilde( 2,:,l,nu);
        bxzeta2 = bxzetatilde(2,:,l,nu);
        for i = 1:Nocc
            fac = VNew(sel(nu),i) * zweight(l)...
                / (zshift(l) - DNew(i));
            W1(:,nu) = W1(:,nu) + ...
                bxzeta1 * bV(:,i) * fac - conj(bzeta1) * bxV(:,i)*conj(fac);
            W2(:,nu) = W2(:,nu) + ...
                bxzeta2 * bV(:,i) * fac - conj(bzeta2) * bxV(:,i)*conj(fac);
        end
    end
    W1(:,nu) = 1/(2i) * W1(:,nu) ;
    W2(:,nu) = 1/(2i) * W2(:,nu) ;
    W(:,nu)  = real((W1(:,nu)) + (W2(:,nu)) + (conj(W2(:,nu))));
    
end

end
