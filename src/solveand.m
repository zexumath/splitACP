function [x1, iter, update] = solveand(GFUNC, x0, b, itertol, alpha)
% Equation to solve is x = g(x) + b;

N = length(x0);
mixdim = 5;
maxite = 50;

rStore = zeros(N,maxite+1);
yStore = zeros(N,maxite+1);
sStore = zeros(N,maxite+1);
sStore(:,1) = x0;


iter = 1;
Fx0 = GFUNC(x0) + b;
rStore(:,iter) = x0 - Fx0;
yStore(:,iter) = rStore(:,iter);
x1 = x0 - alpha(rStore(:,iter));
sStore(:,iter+1) = x1 - x0;
update = norm(sStore(:,iter+1));
iter = iter+1;

while update>itertol && iter < maxite
    x0 = x1;
    Fx0 = GFUNC(x0) + b;
    rStore(:,iter) = x0 - Fx0;
    yStore(:,iter) = rStore(:,iter) - rStore(:,iter-1);
    
    mk = min(iter-1,mixdim);
    Sk = sStore(:,iter:-1:iter-mk+1);
    Yk = yStore(:,iter:-1:iter-mk+1);
    tmpYkrS = Yk\rStore(:,iter);
    sStore(:,iter+1) = -alpha(rStore(:,iter)) + (alpha(Yk)-Sk)*tmpYkrS;
    x1 = x0 + sStore(:,iter+1);
    update = norm(sStore(:,iter+1));
%     update
    iter = iter + 1;
end

end