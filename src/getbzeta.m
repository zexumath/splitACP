function [bzetatilde, bxzetatilde] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatilde,select)

NsCell = opt1D.atom.NsCell;
Npole  = opt1D.Npole;


bzetatilde  = zeros(NsCell,select,Npole,2);
bxzetatilde = zeros(NsCell,select,Npole,2);


for l = 1:Npole
    bzetatilde(:, :, l, 1)     = bV(:,1:Nocc) * zetatilde(1:Nocc,:,l);
    bxzetatilde(:, :, l, 1)    = bxV(:,1:Nocc) * zetatilde(1:Nocc,:,l);
    bzetatilde(:, :, l, 2)     = bV(:,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,:,l);
    bxzetatilde(:, :, l, 2)    = bxV(:,Nocc+1:Ntot) * zetatilde(Nocc+1:Ntot,:,l);
end

end