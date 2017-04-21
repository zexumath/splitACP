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

    
    drhodRb  = zeros(NsGrid,NsCell);
    drhodRbs = zeros(NsGrid,NsCell);
    drhodRbr = zeros(NsGrid,NsCell);
    global tmp0;
    global RELRES;
    RELRES = [];
    tmp0 = zeros(NsGrid,Ntot,NsCell);
    gpfunc = @(x,i)gp(:,i).*x ...
            - b(:,i) * (bx(:,i)'*x)*hs...
            - bx(:,i) * (b(:,i)'*x)*hs;
    for i = 1:NsCell
        gfunc = @(x)gpfunc(x,i);
        drhodRbr(:,i) = adlerwiser1D4m(gfunc,HMatNew,DNew,Vocc,occ, QtProj,opt1D,i);
    end
    
    
%     ita = zeros(NsGrid,Ntot^2);
%     theta = zeros(NsGrid,Ntot^2,Ne);
%     dd  = zeros(Ntot^2,1);
%     DIM = 1;
%     na = Ne;
%     for i = 1:Ntot
%         for j = 1:Ntot
%             ita(:,i+(j-1)*Ntot) = conj(Vocc(:,i)).*(Vocc(:,j));
%             for ia = 1:na
%                 theta(:,i+(j-1)*Ntot,ia) = conj(Vocc(:,j)).*gpfunc(Vocc(:,i),ia);
%             end
%             if(abs(occ(i) - occ(j))>1e-8 ) 
%                 dd(i+(j-1)*Ntot) = (occ(i) - occ(j))/(DNew(i) - DNew(j));
%             else
%                 dd(i+(j-1)*Ntot) = fermidiracderiv(DNew(i),efermi,Tbeta);
%             end
%         end
%     end
%     dd = dd / NsGrid * Ls;
% 
%     drhodRb1 = real(ita * diag(dd) * reshape(sum(theta,1),Ntot^2,na));
    
    eqnsvb = tmp0;
    Npole = 60;
    
    T       = opt1D.temperature;
    Gap     = 0.0;
    DeltaE  = 2.0;
    mu      = efermi;
    
    
    testdPs = zeros(NsGrid,NsGrid,NsCell);
    
    [zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);
    for j = 1:NsCell
        res  = zeros(NsGrid, 1);
        res3 = zeros(NsGrid, NsGrid);
        for l = 1:Npole
            Gl = @(x) VNew *(diag(DNew - zshift(l))\(VNew'*x))* hs;
            tmp1 = zeros(NsGrid, 1);
            tmp3 = zeros(NsGrid, NsGrid);
            for i = 1:Ntot
                tmp2 = gpfunc(VNew(:,i),j);
                tmp1 = tmp1 + Gl( tmp2 ) / (zshift(l) - DNew(i)) .* conj(VNew(:,i));
                tmp3 = tmp3 + Gl( tmp2 ) / (zshift(l) - DNew(i)) * VNew(:,i)';
            end
            res  = res  + zweight(l) * tmp1;
            res3 = res3 + zweight(l) * tmp3;
        end
        res  = 1/(2i)*(res-conj(res));
        res3 = 1/(2i)*(res3-res3');
        drhodRbs(:,j) = res;
        testdPs(:,:,j) = res3;
    end
    drhodRb = drhodRbs + drhodRbr;
    % Then solve equation drhodR = drhodRb + chi_0 v_c drhodR
    
    
    drhodR = zeros(NsGrid,NsCell);
    Cfunc = @(g,index) adlerwiser1D4m(@(x)g.*x, HMatNew, DNew, Vocc,occ, QtProj, opt1D,index);
    Sfunc = @(g)       applychi0s2(opt1D,g, DNew,VNew,efermi,Npole);
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

if(1)
    tic
    
    zeta = eqnsvsc + eqnsvb;
    testdPr = zeros(NsGrid,NsGrid,NsCell);
    for J = 1:NsCell
        for i = 1:Nocc
            testdPr(:,:,J) = testdPr(:,:,J) + occ(i) * VNew(:,i) * zeta(:,i,J)' ...
                + occ(i) * zeta(:,i,J) * VNew(:,i)'; 
        end
    end
end

if(1)
    Vcg = cal_hartree1D(drhodR, opt1D);
    testdPs = real(testdPs + applychi0s3(opt1D,Vcg, DNew,VNew,efermi,Npole));
    
    DFPTtime.constructDiag = toc;
end
%%
if(1)
    testdP = testdPr + testdPs;
    HessDFPT = zeros(NsCell,NsCell);
    for I = 1:NsCell
        for J = 1:NsCell
            HessDFPT(I,J) = HessDFPT(I,J) + sum(gpx(:,I).* drhodR(:,J))*hs;
%             for i = 1:Ne
%                 HessDFPT(I,J) = HessDFPT(I,J)...
%                     + 2*(-1)*hs^2 * sum(VNew(:,i) .* bx(:,I)) * sum(zeta(:,i,J).* b(:,I));
                
%             end
            HessDFPT(I,J) = HessDFPT(I,J)...
                    + 2*(-1)* hs^2 * sum(sum(bx(:,I)' * testdP(:,:,J) * b(:,I)));
        end
        HessDFPT(I,I) = HessDFPT(I,I) + sum(rhoNew.*gpxx(:,I)) * hs;
        for i = 1:Nocc
            HessDFPT(I,I) = HessDFPT(I,I) + 2*(-1) * occ(i)* hs^2 * ( ...
                sum(VNew(:,i) .* bxx(:,I)) * sum(VNew(:,i).*b( :,I)) + ...
                sum(VNew(:,i) .* bx(:, I)).^2 );
        end
    end
    HessDFPT = (HessDFPT + HessDFPT') /2;
    HessDFPT = HessDFPT + HessII;
    DFPTtime.total = toc(DFPTstart);
end
        
            
            


