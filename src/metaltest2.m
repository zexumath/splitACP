% This script implement the FD approach with nlpp term

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
    
    [Energy, force, rho] = scfnlpp4m(opt1D);
end

%%
if(1)
    tic
    forcePerturbed = zeros(NsCell,NsCell);
    rhoPerturbNew = zeros(NsGrid,NsCell);
    EnergyPerturbed = zeros(1,NsCell);
    PPerturb = zeros(NsGrid,NsGrid,NsCell);
    delta = 0.0001;
    for i = 1:NsCell
        opt1Dtmp = opt1D;
        opt1Dtmp.atom.ion(i) = opt1Dtmp.atom.ion(i) + delta;
        opt1Dtmp = initOptnlpp4m(opt1Dtmp);
        fprintf('Atom: \t%d\n',i);
        [EnergyPerturbed(:,i), forcePerturbed(:,i), rhoPerturbNew(:,i)] = scfnlpp4m(opt1Dtmp);
%         [rhoPerturbNew(:,i), ~, ~, tmp, ~, ~, ~] = pptest1(opt1Dtmp);
%         PPerturb(:,:,i) = tmp(:,1:Ne) * tmp(:,1:Ne)';

    end
    
    
    
    
%     HessFD = zeros(NsCell);
%     
%     for i = 1:NsCell
%         for j = 1:NsCell
%             HessFD(i,j) = ...
%                  sum(gpx(:,j).*(rhoPerturbNew(:,i)-rhoNew)/delta ) * hs...
%                 + 2*sum(opt1D.nlppx(:,i).* sum(Vocc,2))...
%                 * sum(opt1D.nlpp(:,i).* sum((VPerturbNew(:,:,j) - Vocc)/delta,2)) * hs^2 ...
%                 + 2*sum(opt1D.nlppx(:,i).* sum((VPerturbNew(:,:,j) - Vocc)/delta,2))...
%                 * sum(opt1D.nlpp(:,i).* sum(Vocc,2)) * hs^2;
%         end
%         diffI = sum(rhoNew.*gpxx(:,i))*hs ...
%             + 2 *hs^2 * sum(sum(Vocc,2) .* opt1D.nlppxx(:,i)) ...
%                 * sum(sum(Vocc,2) .* opt1D.nlpp(:,i)) ...
%             + 2 *hs^2 * sum(sum(Vocc,2) .* opt1D.nlppx(:,i))^2 ;
%         HessFD(i,i) = HessFD(i,i) + diffI;
%     end
%     HessFD = (HessFD + HessFD')/2;
%     eigsFD = sort(eig(HessFD + HessII));
    HessFD = - forcePerturbed + repmat(force,1,NsCell);
    HessFD = HessFD ./delta;
    HessFD = (HessFD + HessFD')/2;
    FDtime = toc;
end