% This script test on if chi_0 matches FD and DFPT with nlpp term
if(1)
%     opt1D = initOptnlpp4m();
    atom = opt1D.atom;
    hs = opt1D.hs;
    Ls = opt1D.Ls;
    NsGrid = opt1D.NsGrid;
    itertol = opt1D.itertol;
    NsCell = atom.NsCell;
    Ne = opt1D.Ne;
    
    rptGrid = (0:NsGrid-1)'*hs;
    
    GptGrid = opt1D.GptGrid;
    
    m0 = opt1D.m0;
    mTotal0 = opt1D.mTotal0;
    mx = opt1D.mx;
    mxx = opt1D.mxx;
    
    gpx = cal_hartree1D(mx,opt1D);
    gpxx = cal_hartree1D(mxx,opt1D);
    
    gp = gpx;
    if(1)
        HessII = zeros(NsCell);
        for i = 1:NsCell
            for j = 1:NsCell
                HessII(i,j) = ...
                    [sum(gpx(:,i).*mx(:,j))]*hs;
                if i==j
                    HessII(i,j) = HessII(i,j)+...
                        [sum(gpxx(:,i).*mTotal0)]*hs;
                end
            end
        end
        HessII = (HessII + HessII') / 2;
    end
    
    [rhoNew, VpotNew, DNew, VNew, occ,efermi, INDNew, HMatNew, iter,update] = metaltest1(opt1D);
    
    id = find(occ>1e-6);
    Nocc = length(id);
    Ntot = length(occ);
    Vocc = VNew(:,1:Nocc);
    
    PcProj = @(x)(Vocc*(Vocc'*x))*hs;
    QcProj = @(x) x - PcProj(x);
    
    PtProj = @(x)(VNew*(VNew'*x))*hs;
    QtProj = @(x) x - PtProj(x);
end


%%
if(1)
    % First compute drhob
    DFPTtime = struct;
    DFPTstart = tic;
    
    b = opt1D.nlpp;
    bx = opt1D.nlppx;
    bxx = opt1D.nlppxx;
    
    temperature = opt1D.temperature;

    Tbeta       = temperature*3.166815d-6;

    
%     drhodRb  = zeros(NsGrid,NsCell);
%     drhodRbs = zeros(NsGrid,NsCell);
    drhodRbr = zeros(NsGrid,NsCell);
    global tmp0;
    global RELRES;
    RELRES = [];
    tmp0 = zeros(NsGrid,Nocc,NsCell);
    gpfunc = @(x,i)gp(:,i).*x ...
            - b(:,i) * (bx(:,i)'*x)*hs...
            - bx(:,i) * (b(:,i)'*x)*hs;
    for i = 1:NsCell
        gfunc = @(x)gpfunc(x,i);
        drhodRbr(:,i) = adlerwiser1D4m(gfunc,HMatNew,DNew,Vocc,occ, QtProj,opt1D,i);
    end
    eqnsvb = tmp0;
    
    ita = zeros(NsGrid,Ntot^2);
    theta = zeros(NsGrid,Ntot^2,Ne);
    dd  = zeros(Ntot^2,1);
    DIM = 1;
    na = Ne;
    for i = 1:Ntot
        for j = 1:Ntot
            ita(:,i+(j-1)*Ntot) = conj(VNew(:,i)).*(VNew(:,j));
            for ia = 1:na
                theta(:,i+(j-1)*Ntot,ia) = conj(VNew(:,j)).*gpfunc(VNew(:,i),ia);
            end
            if(abs(occ(i) - occ(j))>1e-8 ) 
                dd(i+(j-1)*Ntot) = (occ(i) - occ(j))/(DNew(i) - DNew(j));
            else
                dd(i+(j-1)*Ntot) = fermidiracderiv(DNew(i),efermi,Tbeta);
            end
        end
    end
    dd = dd / NsGrid * Ls;

    drhodRbs = real(ita * diag(dd) * reshape(sum(theta,1),Ntot^2,na));
    
    drhodRb = drhodRbs + drhodRbr;
    % Then solve equation drhodR = drhodRb + chi_0 v_c drhodR
    
    drhodR = zeros(NsGrid,NsCell);
    Cfunc = @(g,index) adlerwiser1D4m(@(x)g.*x, HMatNew, DNew, Vocc,occ, QtProj, opt1D,index);
    Sfunc = @(g)       applychi0s1(opt1D,g, DNew,VNew,occ,efermi);
    GFUNC = @(g,index) Cfunc(cal_hartree1D(g, opt1D),index) + Sfunc(cal_hartree1D(g, opt1D));
    
    for j = 1:NsCell
        alpha = 0.5;
        prec = @(x) alpha * x;
        RELRES = [];
        [drhodR(:,j),iter,update] = solveand(@(g)GFUNC(g,j),...
            drhodRb(:,j), drhodRb(:,j), opt1D.pcgtol, prec);
        fprintf('---\n iter = %d, update = %2.15f\n---\n',iter,update);
    end
    eqnsvsc = tmp0;
    
    DFPTtime.response = toc(DFPTstart);
    fprintf('Time consumed for solving equations:\t %10.3g.\n',DFPTtime.response);
end

%%
if(1)
    tic
    Vcg = cal_hartree1D(drhodR, opt1D);
    
    gpV  = zeros(NsGrid,Ntot,NsCell);
    VgpV = zeros(Ntot,Ntot,NsCell);
    for i = 1:Ntot
        for J = 1:NsCell
            gpV(:,i,J) = gpfunc(VNew(:,i),J) + Vcg(:,J) .* VNew(:,i);
            VgpV(:,i,J)  = hs * VNew' * gpV(:,i,J);
        end
    end 
    
    
    
    bV  = hs * b'  * VNew;
    bxV = hs * bx' * VNew;
    
    AB = zeros(NsCell,Ntot^2);
    for i = 1:Ntot
        for a = 1:Ntot
            AB(:,(i-1)*Ntot+a) = bV(:,a).*conj(bxV(:,i));
        end
    end
    dd = dd / hs;
    HessTTs = 2*(-1)*(AB * diag(dd))*reshape(VgpV,Ntot^2,NsCell);
    
    % singular part correct
    
    % regular part
    
    HessTTr = zeros(NsCell,NsCell);
    
    zeta = eqnsvsc + eqnsvb;
    
    for J = 1:NsCell
        HessTTr(:,J) = 2*(-1) * ( ...
         sum( (bV(:,1:Nocc)* diag(occ(1:Nocc))) .* conj(hs*bx' * zeta(:,:,J)),2 ) ...
        +sum( (hs*b'*zeta(:,:,J)).* (conj(bxV(:,1:Nocc)) * diag(occ(1:Nocc))) , 2) );
    end
    
    DFPTtime.constructDiag = toc;
end
%%
if(1)
    HessDFPT = zeros(NsCell,NsCell);
    for I = 1:NsCell
        for J = 1:NsCell
            HessDFPT(I,J) = HessDFPT(I,J) + sum(gpx(:,I).* drhodR(:,J))*hs;
        end
        HessDFPT(I,I) = HessDFPT(I,I) + sum(rhoNew.*gpxx(:,I)) * hs;
        for i = 1:Nocc
            HessDFPT(I,I) = HessDFPT(I,I) + 2*(-1) * occ(i)* hs^2 * ( ...
                sum(VNew(:,i) .* bxx(:,I)) * sum(VNew(:,i).*b( :,I)) + ...
                sum(VNew(:,i) .* bx(:, I)).^2 );
        end
    end
    
    HessDFPT = HessDFPT + HessTTr + HessTTs;
    
    HessDFPT = (HessDFPT + HessDFPT') /2;
    HessDFPT = HessDFPT + HessII;
    DFPTtime.total = toc(DFPTstart);
end
        
            
            


