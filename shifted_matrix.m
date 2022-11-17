function out = shifted_matrix(k,n)
    temp = zeros(n);
    for i=1:n
        if i+k <= n
            temp(i,i+k) = 1;
        end
    end
    idx = 0;
    for i=(n-k+1):n
        idx = idx + 1;
        temp(i, idx) = 1;
    end
    out = temp;
    
end

