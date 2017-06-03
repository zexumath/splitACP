HessTTreg = zeros(NsCell,NsCell);

% ATT    = zeros(NsCell,Nocc);
% BTT    = zeros(Nocc, NsCell);
% CTT    = zeros(NsCell,selectb,NchebNodes);
% DTT    = zeros(selectb,NsCell);
% ETT    = zeros(Nocc,NchebNodes);

cl = DNew(1);
cr = DNew(Nocc);

nodes = cheb_nodes(NchebNodes,[cl,cr]);
[~,ETTr] = lagpoly(DNew(1:Nocc),nodes);

ETTr = diag(occ(1:Nocc))*ETTr;

bxzetabare = zeros(NsCell,selectb,NchebNodes);
bzetabare  = zeros(NsCell,selectb,NchebNodes);

for c = 1:NchebNodes
    bxzetabare(:,:,c) = hs*bx'*eqnresb(:,:,c);
    bzetabare(:,:,c ) = hs*b' *eqnresb(:,:,c);
end


for c = 1:NchebNodes
    AEB1 = bxV(:,1:Nocc)* (repmat(ETTr(:,c),1,NsCell).*bxV(:,1:Nocc)');
    CD1  = (-1)* bzetabare(:,:,c) * b(selb,:);
    AEB2 = bxV(:,1:Nocc)* (repmat(ETTr(:,c),1,NsCell).*bV(:,1:Nocc)');
    CD2  = (-1)* bzetabare(:,:,c) * bx(selb,:);
    AEB3 = bV(:,1:Nocc) * (repmat(ETTr(:,c),1,NsCell).*bxV(:,1:Nocc)');
    CD3  = (-1)* bxzetabare(:,:,c) * b(selb,:);
    AEB4 = bV(:,1:Nocc) * (repmat(ETTr(:,c),1,NsCell).*bV(:,1:Nocc)');
    CD4  = (-1)* bxzetabare(:,:,c) * bx(selb,:);
    
    HessTTreg = HessTTreg + 2*(-1) *(AEB1.*CD1 + AEB2.*CD2 ...
                             + AEB3.*CD3 + AEB4.*CD4 );
end
HessTTreg = (HessTTreg + HessTTreg')/2;

HessTTreg1 = zeros(NsCell,NsCell);

for c = 1:NchebNodes
    AEB1 = bxV(:,1:Nocc) * (diag(ETTr(:,c))*VNew(selb,1:Nocc)');
    AEB2 = bV(:,1:Nocc)  * (diag(ETTr(:,c))*VNew(selb,1:Nocc)');
    HessTTreg1 = HessTTreg1 + 2*(-1) * (AEB1.*bzetabare(:,:,c) + AEB2.*bxzetabare(:,:,c))...
        * gp(selb,:);
end

HessTTreg1 = (HessTTreg1 + HessTTreg1')/2;

% correct so far

Vcg = cal_hartree1D(drhodRACP, opt1D);

bxzeta = zeros(NsCell,select,NchebNodes);
bzeta  = zeros(NsCell,select,NchebNodes);

for c = 1:NchebNodes
    bxzeta(:,:,c) = hs*bx'*eqnres(:,:,c);
    bzeta(:,:,c ) = hs*b' *eqnres(:,:,c);
end

HessTTreg2 = zeros(NsCell,NsCell);

for c = 1:NchebNodes
    AEB1 = bxV(:,1:Nocc) * (diag(ETTr(:,c))*VNew(sel,1:Nocc)');
    AEB2 = bV(:,1:Nocc)  * (diag(ETTr(:,c))*VNew(sel,1:Nocc)');
    HessTTreg2 = HessTTreg2 + 2*(-1) * (AEB1.*bzeta(:,:,c) + AEB2.*bxzeta(:,:,c))...
        * Vcg(sel,:);
end

HessTTreg2 = (HessTTreg2 + HessTTreg2')/2;


%%
% compare with original part

if(0)
    zetatt = zeros(NsGrid,Nocc,NsCell);
    for i = 1:Nocc
        tmp0 = lagpoly(DNew(i),nodes,eqnresb);
        zetatt(:,i,:) = tmp0 *(gp(selb,:).*repmat(Vocc(selb,i),1,NsCell) + (-1)* b(selb,:) * diag(hs*bx'*Vocc(:,i))...
            + (-1)* bx(selb,:) * diag(hs*b'*Vocc(:,i)));
    end
    
    
    testdPACPr = zeros(NsGrid,NsGrid,NsCell);
    for J = 1:NsCell
        for i = 1:Nocc
            testdPACPr(:,:,J) = testdPACPr(:,:,J) + occ(i) * VNew(:,i) * zetatt(:,i,J)' ...
                + occ(i) * zetatt(:,i,J) * VNew(:,i)';
        end
    end
    
    testdPACP = testdPACPr;
    
    
    HessTTACP = zeros(NsCell,NsCell);
    for I = 1:NsCell
        for J = 1:NsCell
            HessTTACP(I,J) = HessTTACP(I,J)...
                + 2*(-1)* hs^2 * sum(sum(bx(:,I)' * (testdPACPr(:,:,J) * b(:,I))));
        end
    end
    HessTTACP = (HessTTACP + HessTTACP') /2;
end


%%

if(1)
    zeta = eqnsvsc + eqnsvb;
    testdPACPr = zeros(NsGrid,NsGrid,NsCell);
    for J = 1:NsCell
        for i = 1:Nocc
            testdPACPr(:,:,J) = testdPACPr(:,:,J) + occ(i) * VNew(:,i) * zeta(:,i,J)' ...
                + occ(i) * zeta(:,i,J) * VNew(:,i)';
        end
    end
end


if(1)
    testdPACP = testdPACPr;
    
    HessACP = zeros(NsCell,NsCell);
    for I = 1:NsCell
        for J = 1:NsCell
            HessACP(I,J) = HessACP(I,J)...
                + 2*(-1)* hs^2 * sum(sum(bx(:,I)' * (testdPACPr(:,:,J) * b(:,I))));
            
        end
    end
    HessACP = (HessACP + HessACP') /2;
end
