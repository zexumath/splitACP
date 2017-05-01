% this script tests on the change of system size

opt1D = initOptnlpp4m();

for Natom = 4:2:12
    opt1D.atom = initAtom1D(Natom);
    opt1D = initOptnlpp4m(opt1D);
    metaltest2;
    metaltest3;
    metaltest4;
    resfilename = sprintf('../res/res0430/metal1Datom%d.mat',atom.NsCell);
    
    eigFD   = eig(HessFD);
    eigDFPT = eig(HessDFPT);
    eigACP  = eig(HessACP);
    
    errorFD  = zeros(opt1D.atom.NsCell,3);
    errorACP = zeros(opt1D.atom.NsCell,3);
    for I = 1:opt1D.atom.NsCell
        errorFD(I,2)    = norm((rhoPerturbNew(:,I) - rho)/delta - drhodR(:,I));
        errorACP(I,1) = norm(drhodRbACP(:,I) - drhodRb(:,I));
        errorACP(I,2) = norm(drhodRACP(:,I)  - drhodR(:,I));
    end
    errorFD(:,3)    = eigFD  - eigDFPT;
    errorACP(:,3)   = eigACP - eigDFPT;
    save(resfilename,'errorFD','errorACP','eigDFPT','eigACP','eigFD','FDtime','DFPTtime','sACPtime');
end
