function res = getdbPb2(opt1D,bV,bxV,I,J,DNew,VNew,Nocc,zshift,zweight,bzetatilde, bxzetatilde,sel)
Npole   = opt1D.Npole;
hs      = opt1D.hs;
select  = length(sel);
NsGrid  = opt1D.NsGrid;
Ntot    = length(DNew);
res     = 0;

b       = opt1D.nlpp;
bx      = opt1D.nlppx;


W1 = zeros(1,select);
W2 = zeros(1,select);
W  = zeros(1,select);
for nu = 1:select
    if b(sel(nu),J)==0 && bx(sel(nu),J) == 0
        continue
    else
        for l = 1:Npole
            bzeta1  = bzetatilde( 1,I,l,nu);
            bxzeta1 = bxzetatilde(1,I,l,nu);
            bzeta2  = bzetatilde( 2,I,l,nu);
            bxzeta2 = bxzetatilde(2,I,l,nu);
            for i = 1:Nocc
                %             fac = (VNew(sel(nu),i)*gp(sel(nu),J) ...
                %                 +(-1)*b(sel(nu),J)*bxV(J,i) + (-1)*bx(sel(nu),J)*bV(J,i))...
                %                 * zweight(l)/ (zshift(l) - DNew(i));
                %
                
                fac = ((-1)*b(sel(nu),J)*bxV(J,i) + (-1)*bx(sel(nu),J)*bV(J,i))...
                    * zweight(l)/ (zshift(l) - DNew(i));
                
                W1(:,nu) = W1(:,nu) + ...
                    bxzeta1 * bV(I,i) * fac - conj(bzeta1) * bxV(I,i)*conj(fac);
                W2(:,nu) = W2(:,nu) + ...
                    bxzeta2 * bV(I,i) * fac - conj(bzeta2) * bxV(I,i)*conj(fac);
            end
        end
        W1(:,nu) = 1/(2i) * W1(:,nu);
        W2(:,nu) = 1/(2i) * W2(:,nu);
        res = res + real((W1(:,nu)) + (W2(:,nu)) + (conj(W2(:,nu))));
    end
end

end
