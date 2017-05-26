function specplot(X,sigma,range,linspec,linewidth)
if nargin <=1
    sigma =0.001;
end
if nargin <=2
    range = [-10*sigma,X(end)+10*sigma];
end
if nargin <=3
    linspec = 'b-';
end
if nargin <=4
    linewidth = 1;
end


gaussian1D = @(x,mu,sigma) 1/sqrt(2*pi)/sigma * exp(-(x-mu).^2/2/sigma^2);
func = @(x) sum(gaussian1D(x,X,sigma))/length(X);
% gaussian1D(0,X,sigma)
% fplot((@(x)gaussian1D(x,X(1),sigma)),range);

figure(1)
xgrid = (0:0.001:1) * (range(2) - range(1)) + range(1);
y = zeros(1001,1);
for i = 1:1001
    y(i) = func(xgrid(i));
end
plot(xgrid,y,linspec,'LineWidth',linewidth);

