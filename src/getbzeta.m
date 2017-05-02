function [bzetatilde, bxzetatilde] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatilde,select)

NsCell = opt1D.atom.NsCell;
Npole  = opt1D.Npole;


bzetatilde  = zeros(2,NsCell,Npole,select);
bxzetatilde = zeros(2,NsCell,Npole,select);


for nu = 1:select
    for l = 1:Npole
        for I = 1:NsCell
            bzetatilde(1, I, l, nu)     = bV(I,1:Nocc) * zetatilde(1:Nocc,l,nu);
            bxzetatilde(1, I, l, nu)    = bxV(I,1:Nocc) * zetatilde(1:Nocc,l,nu);
            bzetatilde(2, I, l, nu)     = bV(I,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,l,nu);
            bxzetatilde(2, I, l, nu)    = bxV(I,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,l,nu);
        end
    end
end

end