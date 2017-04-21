% test on change of extra orbitals ntot

opt1D = initOptnlpp4m();
metaltest3;
for testratio = 0.9-0.1:0.2
    opt1D.NeExtra = ceil(opt1D.Ne * testratio);
    metaltest4;
    
    eigDFPT = eig(HessDFPT);
    eigACP  = eig(HessACP);
    
    errorACP = zeros(opt1D.atom.NsCell,3);
    for I = 1:opt1D.atom.NsCell
        errorACP(I,1) = norm(drhodRbACP(:,I) - drhodRb(:,I));
        errorACP(I,2) = norm(drhodRACP(:,I)  - drhodR(:,I));
    end
    errorACP(:,3)   = eigACP - eigDFPT;
    
    resfilename = sprintf('./result041917/metal1Dratio%1.1s.mat',testratio);
    save(resfilename,'errorACP','eigDFPT','eigACP','DFPTtime','sACPtime');
end