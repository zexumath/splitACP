function [result, weight] = lagpoly(x,nodes,values)

if(nargin ==2)
    result = 0;
    Nc = length(nodes);
    weight = zeros(length(x),Nc);
    for i=1:Nc
        tmp = 1;
        if(1)
            for j=1:Nc
                if j~=i
                    tmp = tmp .* (x - nodes(j)) / (nodes(i)-nodes(j));
                end
            end
            weight(:,i) = tmp;
        end
    end
else
    result = 0;
    Nc = length(nodes);
    weight = zeros(Nc,1);
    for i=1:Nc
        tmp = 1;
        if(1)
            for j=1:Nc
                if j~=i
                    tmp = tmp * (x - nodes(j)) / (nodes(i)-nodes(j));
                end
            end
            weight(i) = tmp;
            result = result + tmp * values(:,:,i);
        end
    end
end
end