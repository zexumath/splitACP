

HessTTsing = zeros(NsCell,NsCell);

ETT = 1./(-repmat(DNew(1:Nocc),1,Npole) + repmat(conj(zshift'),Nocc,1));

ETT = ETT*diag(zweight);
[bzetatildebare, bxzetatildebare] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatildeb,selectb);

for l = 1:Npole
    
    %(3.13) line three
    
    AEB1 = bxV(:,1:Nocc)* (repmat(ETT(:,l),1,NsCell).*bxV(:,1:Nocc)');
    CD1  = (-1)* bzetatildebare(:,:,l,1) * b(selb,:);
    AEB2 = bxV(:,1:Nocc)* (repmat(ETT(:,l),1,NsCell).*bV(:,1:Nocc)');
    CD2  = (-1)* bzetatildebare(:,:,l,1) * bx(selb,:);
    
    AEB3 = bV(:,1:Nocc) * (repmat(-conj(ETT(:,l)),1,NsCell).*bxV(:,1:Nocc)');
    CD3  = (-1)* conj(bxzetatildebare(:,:,l,1)) * b(selb,:);
    AEB4 = bV(:,1:Nocc) * (repmat(-conj(ETT(:,l)),1,NsCell).*bV(:,1:Nocc)');
    CD4  = (-1)* conj(bxzetatildebare(:,:,l,1)) * bx(selb,:);
    
    HessTTsing = HessTTsing + 1/(2i)* 2*(-1)* ...
        (AEB1.*CD1 + AEB2.*CD2 + AEB3.*CD3 + AEB4.*CD4 );
    
    %(3.13) line two with h.c. included
    
    AEB1 = bxV(:,1:Nocc)* (repmat(ETT(:,l),1,NsCell).*bxV(:,1:Nocc)');
    CD1  = (-1)* bzetatildebare(:,:,l,2) * b(selb,:);
    AEB2 = bxV(:,1:Nocc)* (repmat(ETT(:,l),1,NsCell).*bV(:,1:Nocc)');
    CD2  = (-1)* bzetatildebare(:,:,l,2) * bx(selb,:);
    
    AEB3 = bV(:,1:Nocc) * (repmat(-conj(ETT(:,l)),1,NsCell).*bxV(:,1:Nocc)');
    CD3  = (-1)* conj(bxzetatildebare(:,:,l,2)) * b(selb,:);
    AEB4 = bV(:,1:Nocc) * (repmat(-conj(ETT(:,l)),1,NsCell).*bV(:,1:Nocc)');
    CD4  = (-1)* conj(bxzetatildebare(:,:,l,2)) * bx(selb,:);
    
    HessTTsing = HessTTsing +  2*(-1)* ...
        imag((AEB1.*CD1 + AEB2.*CD2 + AEB3.*CD3 + AEB4.*CD4 ) );
    
    % correct
end
HessTTsing = real(HessTTsing + HessTTsing')/2;

%%

if(1)
    tic
    %singular part dealt seperately.
    testdPACP = testdPACPr;
    
    [bzetatilde, bxzetatilde] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatilde,select);
    [bzetatildebare, bxzetatildebare] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatildeb,selectb);
    HessACP = zeros(NsCell,NsCell);
    for I = 1:NsCell
        for J = 1:NsCell
            HessACP(I,J) = HessACP(I,J)...
                + 2*(-1) * getdbPb2(opt1D,bV,bxV,I,J,DNew,VNew,Nocc,zshift,zweight,bzetatildebare, bxzetatildebare,selb);            
        end

    end
    HessACP = (HessACP + HessACP') /2;
end