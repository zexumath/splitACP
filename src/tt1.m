
metaltest2;
eigFD = eig(HessFD);


metaltest3;
eigDFPT = eig(HessDFPT);


metaltest4;
eigACP = eig(HessACP);

error = [eigFD - eigDFPT,eigACP - eigDFPT]

