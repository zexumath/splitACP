function nodes = cheb_nodes(n,range)
    % return chebyshev nodes 
    nodes = zeros(1,n);
    for i=1:n
        nodes(i) = -cos((2.0*i-1.0)*pi/(2.0*n));
    end
    nodes = range(2) + (1.0 + nodes)*(range(1)-range(2))/2;
end