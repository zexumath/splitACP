function result = lagpoly(x,nodes,values)
    result = 0;
    N = length(nodes);
    for i=1:N
        tmp = 1;
        if(1)
            for j=1:N
                if j~=i
                    tmp = tmp * (x - nodes(j)) / (nodes(i)-nodes(j));
                end
            end
            result = result + tmp * values(:,:,i);
        end
    end
end