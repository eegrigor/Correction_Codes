function h = base_to_parity_check(B, z)
    [m,n] = size(B);
    H = zeros(z, z*n);
    temp = zeros(z,z);

    for i=1:m
        temp_row = zeros(z,z);
        for j=1:n
            k = B(i,j);
            if k == -1
                temp = zeros(z,z);
                temp_row = cat(2,temp_row, temp);
            else
                temp = shifted_matrix(k, z);
                temp_row = cat(2,temp_row, temp);
            end


        end
        temp_row = temp_row(:, z+1:end);
        H = cat(1,H, temp_row);
    end
    H = H(z+1:end, :);
    h=H;
end

