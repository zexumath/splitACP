% This script test on the error of split ACP
% Natom = 60; temperature = 2000;
% Npole = 20; NchebNodes = 4;

% opt1D = initOptnlpp4m();
% 
% metaltest2;
% metaltest6;
% 

% errACPspSV = zeros(NsCell,2,6);

idxxx = 1;
for testratio = 0.9:-0.1:0.5
    opt1D.NeExtra = ceil(opt1D.Ne * testratio);
    metaltest5;
    geterror;
    
    idxxx = idxxx + 1;
    errACPspSV(:,:,idxxx) = errorACPsp(:,2:3);
end
