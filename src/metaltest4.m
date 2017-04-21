% This script test on if chi_0 matches FD and ACP with nlpp term
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
    
    %% compute drhodRbCSr
    [RHSsel,sel, drhodRbACPr, CHI0gCStmp, W, PI,sACPtime,pivot,indrand,rphase,eqnres,dpsib] = ...
        compresschi01Dtest(opt1D,gpfunc,Vocc,occ,DNew,HMatNew,QtProj,NchebNodes,sACPtime);
    
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
    gpfunc1 = @(x,i,mu)gp(mu,i)*x(mu) ...
            - b(mu,i) * (bx(:,i)'*x)*hs...
            - bx(mu,i) * (b(:,i)'*x)*hs;
    
    [drhodRbACPs] = generatew2(opt1D,VNew,DNew,Nocc,zshift,zweight,gpfunc1,RHSsel,sel);
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
    fprintf('Time consumed for solving equations:\t %10.3g.\n',sACPtime.response);
end


if(1)
    
    bV  = zeros(NsCell,Ntot);
    bxV = zeros(NsCell,Ntot);
    
    for I = 1:NsCell
        for i = 1:Ntot
            bV(I,i)     = sum(b(:,I).*VNew(:,i))*hs;
            bxV(I,i)    = sum(bx(:,I).*VNew(:,i))*hs;
        end
    end
end



if(1)
    tic
    
    zeta = eqnsvsc + eqnsvb;
    testdPACPr = zeros(NsGrid,NsGrid,NsCell);
    for J = 1:NsCell
        for i = 1:Nocc
            testdPACPr(:,:,J) = testdPACPr(:,:,J) + occ(i) * VNew(:,i) * zeta(:,i,J)' ...
                + occ(i) * zeta(:,i,J) * VNew(:,i)'; 
        end
    end
    
    Vcg = cal_hartree1D(drhodRACP, opt1D);
%     [Ws, Wmat] = generatew(opt1D,VNew,DNew,zshift,zweight,RHSsel,sel);
%     for J = 1:NsCell
%         for nu = 1:select
%             testdPACPs(:,:,J) = testdPACPs(:,:,J) + Wmat(:,:,nu) * Vcg(sel(nu),J);
%         end
%     end
    
    
    
    sACPtime.constructDiag = toc;
end


if(1)
    tic
    
    testdPACP = testdPACPr; % singular part dealt directly
    HessACP = zeros(NsCell,NsCell);
    for I = 1:NsCell
        for J = 1:NsCell
            HessACP(I,J) = HessACP(I,J) + sum(gpx(:,I).* drhodRACP(:,J))*hs;
%             for i = 1:Ne
%                 HessDFPT(I,J) = HessDFPT(I,J)...
%                     + 2*(-1)*hs^2 * sum(VNew(:,i) .* bx(:,I)) * sum(zeta(:,i,J).* b(:,I));
                
%             end
            HessACP(I,J) = HessACP(I,J)...
                    + 2*(-1)* hs^2 * sum(sum(bx(:,I)' * testdPACP(:,:,J) * b(:,I)));
            HessACP(I,J) = HessACP(I,J)...
                    + (-1) * getdbPb(opt1D,bV,bxV,I,J,DNew,VNew,zshift,zweight,Vcg,gpfunc);
        end
        HessACP(I,I) = HessACP(I,I) + sum(rhoNew.*gpxx(:,I)) * hs;
        for i = 1:Nocc
            HessACP(I,I) = HessACP(I,I) + 2*(-1) * occ(i)* hs^2 * ( ...
                sum(VNew(:,i) .* bxx(:,I)) * sum(VNew(:,i).*b( :,I)) + ...
                sum(VNew(:,i) .* bx(:, I)).^2 );
        end
    end
    HessACP = (HessACP + HessACP') /2;
    HessACP = HessACP + HessII;
    
    sACPtime.integral = toc;
    sACPtime.total = toc(ACPstart);
end

        
            
            


