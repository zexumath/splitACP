% This script test on if chi_0 matches FD and ACP with nlpp term
if(1)
    opt1D = initOptnlpp4m();
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
    Docc = DNew(1:Nocc);
    
    %     PcProj = @(x)(Vocc*(Vocc'*x))*hs;
    %     QcProj = @(x) x - PcProj(x);
    
    PtProj = @(x)(VNew*(VNew'*x))*hs;
    QtProj = @(x) x - PtProj(x);
    
end


%%
if(1)
    % First compute drhob
    sACPtime = struct;
    sACPtime.colsel = 0;
    sACPtime.soleqn = 0;
    sACPtime.reconstruct = 0;
    %     CStime.proj = 0;
    sACPtime.iter = [];
    sACPtime.total = 0;
    
    ACPstart = tic;
    
    b = opt1D.nlpp;
    bx = opt1D.nlppx;
    bxx = opt1D.nlppxx;
    
    temperature = opt1D.temperature;
    NchebNodes  = opt1D.NchebNodes;
    Tbeta       = temperature*3.166815d-6;
    
    
    %     drhodRb = zeros(NsGrid,NsCell);
    %     drhodRb1 = zeros(NsGrid,NsCell);
    %     drhodRb2 = zeros(NsGrid,NsCell);
    
    tmp0 = zeros(NsGrid,Ntot,NsCell);
    
    gpfunc = @(x,i)gp(:,i).*x ...
        - b(:,i) * (bx(:,i)'*x)*hs...
        - bx(:,i) * (b(:,i)'*x)*hs;
    
    if(1)
        
        bV  = zeros(NsCell,Ntot);
        bxV = zeros(NsCell,Ntot);
        
%         for I = 1:NsCell
%             for i = 1:Ntot
%                 bV(I,i)     = sum(b(:,I).*VNew(:,i))*hs;
%                 bxV(I,i)    = sum(bx(:,I).*VNew(:,i))*hs;
%             end
%         end
        bV = b'*VNew*hs;
        bxV = bx'*VNew*hs;
    end
    
    %% compute drhodRbCSr
    [RHSselb,selb, drhodRbACPr, CHI0gCStmp, W, PI,sACPtime,pivot,indrand,rphase,eqnresb,dpsib] = ...
        compresschi01Dtest(opt1D,gpfunc,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,sACPtime);
    selectb = length(selb);
    sACPtime.regular = toc(ACPstart);
    
    eqnsvb = dpsib;
    
    %% compute drhodRbCSs
    
    singularstart = tic;
    
    Npole   = opt1D.Npole;
    T       = opt1D.temperature;
    Gap     = 0.0;
    DeltaE  = max(2.0, DNew(end) - DNew(1));
    mu      = efermi;
    
    
    
    [zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);
    %     drhodRbACPs = zeros(NsGrid,NsCell);
    %     for j = 1:NsCell
    %         res = zeros(NsGrid, 1);
    %         for l = 1:Npole
    %             Gl = @(x) VNew *((DNew - zshift(l)).\(VNew'*x))* hs;
    %             tmp1 = zeros(NsGrid, 1);
    %             for i = 1:Ntot
    %                 tmp2 = gpfunc(VNew(:,i),j);
    %                 tmp1 = tmp1 + Gl( tmp2 ) / (zshift(l) - DNew(i)) .* conj(VNew(:,i));
    %             end
    %             res = res + zweight(l) * tmp1;
    %         end
    %         res = 1/(2i)*(res-conj(res));
    %         drhodRbACPs(:,j) = res;
    %     end
    %
    % test for using ISDF
    
    zetatildeb = zeros(Ntot,selectb,Npole);
    
    for l = 1:Npole
        zetatildeb(:,:,l) = (diag(DNew - zshift(l))\ (VNew'*RHSselb* hs));
    end
    
    
    W0 = generatew3(opt1D,VNew,DNew,Nocc,zshift,zweight,RHSselb,selb);
    [drhodRbACPs] = generatew2(opt1D,bV,bxV,VNew,DNew,Nocc,zshift,zweight,zetatildeb,selb);
    
    drhodRbACPs = drhodRbACPs + W0*gp(selb,:);
    %
    
    sACPtime.singular = toc(singularstart);
    
    drhodRbACP = drhodRbACPs + drhodRbACPr;
    
    sACPtime.iter(1) = toc(ACPstart);
    
    %% Then solve equation drhodR = drhodRb + chi_0 v_c drhodR
    Npole = opt1D.Npole;
    
    maxiterACP = 5;
    drhodRACPold = drhodRbACP;
    Y0 = cal_hartree1D(drhodRbACP,opt1D);
    
    for iter = 1:maxiterACP
        tmpstart = tic;
        Y = cal_hartree1D(drhodRACPold,opt1D);
        if(iter>1)
            [sel, Wnew, PInew,sACPtime,pivot,indrand,rphase,eqnres,dpsisc,RHSsel] = ...
                compresschi01Dr(opt1D,Y,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,...
                sACPtime,pivot,indrand,rphase,select,eqnres);
        else
            [sel, Wnew, PInew,sACPtime,pivot,indrand,rphase,eqnres,dpsisc,RHSsel] = ...
                compresschi01Dr(opt1D,Y,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,sACPtime);
        end
        select = length(sel);
        
        Ws = generatew3(opt1D,VNew,DNew,Nocc,zshift,zweight,RHSsel,sel);
        Wtot = Wnew + Ws;
        drhodRACP = drhodRbACP + Wtot * ...
            ( (eye(length(sel)) - PInew(cal_hartree1D(Wtot,opt1D)))\...
            PInew(Y0));
        
        errnrm = norm(drhodRACP - drhodRACPold)./norm(drhodRACPold);
        fprintf('iter = %d,\t||CHIgCSnew1 - CHIgCSnew||_2\t= %3.15f \n',iter,errnrm);
        if errnrm < opt1D.requiredeps
            break;
        end
        drhodRACPold = drhodRACP;
        sACPtime.iter(iter+1) = toc(tmpstart);
    end
    eqnsvsc = dpsisc;
    
    sACPtime.response = toc(ACPstart);
    fprintf('Time consumed for solving Dyson equations:\t %10.3g.\n',sACPtime.response);
end

%%

if(1)
    tic
    
    zetatilde = zeros(Ntot,select,Npole);
    
    for l = 1:Npole
        zetatilde(:,:,l) = (diag(DNew - zshift(l))\ (VNew'*RHSsel* hs));
    end

    % first deal with regular part.
    
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
    
    % then singular part
    HessTTsing = zeros(NsCell,NsCell);

    ETTs = 1./(-repmat(DNew(1:Nocc),1,Npole) + repmat(conj(zshift'),Nocc,1));
    
    ETTs = ETTs*diag(zweight);
    [bzetatildebare, bxzetatildebare] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatildeb,selectb);
    
    for l = 1:Npole
        
%         (3.13) line three
        
        AEB1 = bxV(:,1:Nocc)* (repmat(ETTs(:,l),1,NsCell).*bxV(:,1:Nocc)');
        CD1  = (-1)* bzetatildebare(:,:,l,1) * b(selb,:);
        AEB2 = bxV(:,1:Nocc)* (repmat(ETTs(:,l),1,NsCell).*bV(:,1:Nocc)');
        CD2  = (-1)* bzetatildebare(:,:,l,1) * bx(selb,:);
        
        AEB3 = bV(:,1:Nocc) * (repmat(-conj(ETTs(:,l)),1,NsCell).*bxV(:,1:Nocc)');
        CD3  = (-1)* conj(bxzetatildebare(:,:,l,1)) * b(selb,:);
        AEB4 = bV(:,1:Nocc) * (repmat(-conj(ETTs(:,l)),1,NsCell).*bV(:,1:Nocc)');
        CD4  = (-1)* conj(bxzetatildebare(:,:,l,1)) * bx(selb,:);
        
        HessTTsing = HessTTsing + 1/(2i)* 2*(-1)* ...
            (AEB1.*CD1 + AEB2.*CD2 + AEB3.*CD3 + AEB4.*CD4 );
        
%         (3.13) line two with h.c. included
        
        AEB1 = bxV(:,1:Nocc)* (repmat(ETTs(:,l),1,NsCell).*bxV(:,1:Nocc)');
        CD1  = (-1)* bzetatildebare(:,:,l,2) * b(selb,:);
        AEB2 = bxV(:,1:Nocc)* (repmat(ETTs(:,l),1,NsCell).*bV(:,1:Nocc)');
        CD2  = (-1)* bzetatildebare(:,:,l,2) * bx(selb,:);
        
        AEB3 = bV(:,1:Nocc) * (repmat(-conj(ETTs(:,l)),1,NsCell).*bxV(:,1:Nocc)');
        CD3  = (-1)* conj(bxzetatildebare(:,:,l,2)) * b(selb,:);
        AEB4 = bV(:,1:Nocc) * (repmat(-conj(ETTs(:,l)),1,NsCell).*bV(:,1:Nocc)');
        CD4  = (-1)* conj(bxzetatildebare(:,:,l,2)) * bx(selb,:);
        
        HessTTsing = HessTTsing +  2*(-1)* ...
            imag((AEB1.*CD1 + AEB2.*CD2 + AEB3.*CD3 + AEB4.*CD4 ) );
        
%         correct
    end
    HessTTsing = real(HessTTsing + HessTTsing')/2;
    
    
    [bzetatilde, bxzetatilde] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatilde,select);
    [bzetatildeb, bxzetatildeb] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatildeb,selectb);
    
    HessEasy = zeros(NsCell,NsCell);
    for I = 1:NsCell
        WVcg = getdbPb1(opt1D,bV,bxV,I,DNew,VNew,Nocc,zshift,zweight,bzetatilde, bxzetatilde,sel);
        Wb   = getdbPb1(opt1D,bV,bxV,I,DNew,VNew,Nocc,zshift,zweight,bzetatildeb, bxzetatildeb,selb);
        for J = 1:NsCell
            HessEasy(I,J) = HessEasy(I,J) + sum(gpx(:,I).* drhodRACP(:,J))*hs;
            HessEasy(I,J) = HessEasy(I,J)...
                + 2*(-1) * WVcg * Vcg(sel,J)...
                + 2*(-1) * Wb   * gp(selb,J) ;
        end
        HessEasy(I,I) = HessEasy(I,I) + sum(rhoNew.*gpxx(:,I)) * hs;
        for i = 1:Nocc
            HessEasy(I,I) = HessEasy(I,I) + 2*(-1) * occ(i)* hs^2 * ( ...
                sum(VNew(:,i) .* bxx(:,I)) * sum(VNew(:,i).*b( :,I)) + ...
                sum(VNew(:,i) .* bx(:, I)).^2 );
        end
    end
    HessEasy = (HessEasy + HessEasy') /2;
    HessACP = HessEasy + HessTTreg + HessTTreg1 + HessTTreg2 + HessTTsing + HessII;
    
    sACPtime.integral = toc;
    sACPtime.total = toc(ACPstart);
end






