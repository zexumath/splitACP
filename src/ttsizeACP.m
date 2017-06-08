% This script test on the scale of split ACP and ACP. 

opt1D = initOptnlpp4m();

for Natom = 10:10:100
    opt1D.atom = initAtom1D(Natom);
    opt1D = initOptnlpp4m(opt1D);
    metaltest2;
    metaltest5;
    metaltest6;
    metaltest7;
    
    eigFD   = eig(HessFD);
    eigDFPT = eig(HessDFPT);
    eigACP  = eig(HessACP);
    eigACPsp= eig(HessACPsp);
    
    errorFD  = zeros(opt1D.atom.NsCell,3);
    errorACP = zeros(opt1D.atom.NsCell,3);
    errorACPsp = zeros(opt1D.atom.NsCell,3);
    for I = 1:opt1D.atom.NsCell
        errorFD(I,2)  = norm((rhoPerturbNew(:,I) - rho)/delta - drhodR(:,I));
        errorACP(I,1) = norm(drhodRbACP(:,I) - drhodRb(:,I));
        errorACP(I,2) = norm(drhodRACP(:,I)  - drhodR(:,I));
        errorACPsp(I,1) = norm(drhodRbACPsp(:,I) - drhodRb(:,I));
        errorACPsp(I,2) = norm(drhodRACPsp(:,I)  - drhodR(:,I));
    end
    errorFD(:,3)    = eigFD  - eigDFPT;
    errorACP(:,3)   = eigACP - eigDFPT;
    errorACPsp(:,3) = eigACPsp - eigDFPT;
    resfilename = sprintf('../res/res060717/metal1DT0%d.mat',atom.NsCell);
    save(resfilename,'errorFD','errorACP','errorACP',...
        'eigDFPT','eigACP','eigFD','eigACPsp',...
        'FDtime','DFPTtime','spACPtime','ACPtime');
end
