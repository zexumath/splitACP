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
    spACPtime = struct;
    spACPtime.colsel = 0;
    spACPtime.soleqn = 0;
    spACPtime.reconstruct = 0;
    %     CStime.proj = 0;
    spACPtime.iter = [];
    spACPtime.total = 0;
    
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
    [RHSselb,selb, drhodRbACPr, CHI0gCStmp, W, PI,spACPtime,pivot,indrand,rphase,eqnresb,dpsib] = ...
        compresschi01Dtest(opt1D,gpfunc,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,spACPtime);
    selectb = length(selb);
    spACPtime.regular = toc(ACPstart);
    
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
    
    spACPtime.singular = toc(singularstart);
    
    drhodRbACP = drhodRbACPs + drhodRbACPr;
    
    spACPtime.iter(1) = toc(ACPstart);
    
    %% Then solve equation drhodR = drhodRb + chi_0 v_c drhodR
    Npole = opt1D.Npole;
    
    maxiterACP = 5;
    drhodRACPold = drhodRbACP;
    Y0 = cal_hartree1D(drhodRbACP,opt1D);
    
    for iter = 1:maxiterACP
        tmpstart = tic;
        Y = cal_hartree1D(drhodRACPold,opt1D);
        if(iter>1)
            [sel, Wnew, PInew,spACPtime,pivot,indrand,rphase,eqnres,dpsisc,RHSsel] = ...
                compresschi01Dr(opt1D,Y,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,...
                spACPtime,pivot,indrand,rphase,select,eqnres);
        else
            [sel, Wnew, PInew,spACPtime,pivot,indrand,rphase,eqnres,dpsisc,RHSsel] = ...
                compresschi01Dr(opt1D,Y,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,spACPtime);
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
        spACPtime.iter(iter+1) = toc(tmpstart);
    end
    eqnsvsc = dpsisc;
    
    spACPtime.response = toc(ACPstart);
    fprintf('Time consumed for solving Dyson equations:\t %10.3g.\n',spACPtime.response);
end

%%

if(1)
    constructDiagStart = tic;
    
    zetatilde = zeros(Ntot,select,Npole);
    
    for l = 1:Npole
        zetatilde(:,:,l) = (diag(DNew - zshift(l))\ (VNew'*RHSsel* hs));
    end
    
    
    zeta = eqnsvsc + eqnsvb;
    testdPACPr = zeros(NsGrid,NsGrid,NsCell);
    for J = 1:NsCell
        for i = 1:Nocc
            testdPACPr(:,:,J) = testdPACPr(:,:,J) + occ(i) * VNew(:,i) * zeta(:,i,J)' ...
                + occ(i) * zeta(:,i,J) * VNew(:,i)';
        end
    end
    
%     testdPACPs1 = zeros(NsGrid,NsGrid,NsCell);
%     
    Vcg = cal_hartree1D(drhodRACP, opt1D);
%     [~, Wmatb] = generatew3(opt1D,VNew,DNew,Nocc,zshift,zweight,RHSselb,selb);
%     [~, Wmat ] = generatew3(opt1D,VNew,DNew,Nocc,zshift,zweight,RHSsel ,sel );
%     
%     tic
%     for J = 1:NsCell
%         for nu = 1:selectb
%             testdPACPs1(:,:,J) = testdPACPs1(:,:,J) ...
%                 + Wmatb(:,:,nu) * gp(selb(nu),J);
%             
%             if(b(selb(nu),J)~=0)
%                 dPdRACPb1  = zeros(NsGrid,NsGrid);
%                 dPdRACPb2  = zeros(NsGrid,NsGrid);
%                 
%                 
%                 for l = 1:Npole
%                     zeta1 = VNew(:,1:Nocc) * ((DNew(1:Nocc) - zshift(l)).\ (VNew(:,1:Nocc)'*RHSselb(:,nu)* hs));
%                     zeta2 = VNew(:,Nocc+1:Ntot) * ((DNew(Nocc+1:Ntot) - zshift(l)).\ (VNew(:,Nocc+1:Ntot)'*RHSselb(:,nu)* hs));
%                     for i = 1:Nocc
%                         fac = (-1)*(b(selb(nu),J)*bxV(J,i) + bx(selb(nu),J)*bV(J,i)) * zweight(l)...
%                             / (zshift(l) - DNew(i));
%                         dPdRACPb1 = dPdRACPb1 + ...
%                             ( zeta1*VNew(:,i)') * fac;
%                         dPdRACPb2 = dPdRACPb2 + ...
%                             ( zeta2*VNew(:,i)') * fac;
%                     end
%                 end
%                 
%                 dPdRACPb1 = 1/(2i)* (dPdRACPb1 - dPdRACPb1');
%                 dPdRACPb2 = 1/(2i)* (dPdRACPb2 - dPdRACPb2');
%                 
%                 testdPACPs1(:,:,J) = testdPACPs1(:,:,J)+...
%                     real(dPdRACPb1 + dPdRACPb2 + (dPdRACPb2)');
%             end
%         end
%     end
%     sACPtime.constructDiagNLPP = toc;
%     
%     testdPACPs2 = zeros(NsGrid,NsGrid,NsCell);
% 
%     
%     for J = 1:NsCell
%         for nu = 1:select
%             testdPACPs2(:,:,J) = testdPACPs2(:,:,J) ...
%                 + Wmat(:,:,nu) * Vcg(sel(nu),J);
%         end
%     end
%     
%     
%     
    spACPtime.constructDiag = toc(constructDiagStart);
end


if(1)
    tic
    %singular part dealt seperately.
    testdPACP = testdPACPr;
    
    [bzetatilde, bxzetatilde] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatilde,select);
    [bzetatildebare, bxzetatildebare] = getbzeta(opt1D,bV,bxV,Nocc,Ntot,zetatildeb,selectb);
    
    
    HessspACP = zeros(NsCell,NsCell);
    for I = 1:NsCell
        WVcg = getdbPb1(opt1D,bV,bxV,I,DNew,VNew,Nocc,zshift,zweight,bzetatilde, bxzetatilde,sel);
        Wb   = getdbPb1(opt1D,bV,bxV,I,DNew,VNew,Nocc,zshift,zweight,bzetatildebare,bxzetatildebare,selb);
        for J = 1:NsCell
            HessspACP(I,J) = HessspACP(I,J) + sum(gpx(:,I).* drhodRACP(:,J))*hs;
            %             for i = 1:Ne
            %                 HessDFPT(I,J) = HessDFPT(I,J)...
            %                     + 2*(-1)*hs^2 * sum(VNew(:,i) .* bx(:,I)) * sum(zeta(:,i,J).* b(:,I));
            
            %             end
            HessspACP(I,J) = HessspACP(I,J)...
                + 2*(-1)* hs^2 * sum(sum(bx(:,I)' * (testdPACPr(:,:,J) * b(:,I))));
            HessspACP(I,J) = HessspACP(I,J)...
                + 2*(-1) * WVcg * Vcg(sel,J)...
                + 2*(-1) * Wb   * gp(selb,J) ...
                + 2*(-1) * getdbPb2(opt1D,bV,bxV,I,J,DNew,VNew,Nocc,zshift,zweight,bzetatildebare, bxzetatildebare,selb);
            
            
            
        end
        HessspACP(I,I) = HessspACP(I,I) + sum(rhoNew.*gpxx(:,I)) * hs;
        for i = 1:Nocc
            HessspACP(I,I) = HessspACP(I,I) + 2*(-1) * occ(i)* hs^2 * ( ...
                sum(VNew(:,i) .* bxx(:,I)) * sum(VNew(:,i).*b( :,I)) + ...
                sum(VNew(:,i) .* bx(:, I)).^2 );
        end
    end
    HessspACP = (HessspACP + HessspACP') /2;
    HessspACP = HessspACP + HessII;
    
    spACPtime.integral = toc;
    spACPtime.total = toc(ACPstart);
end






