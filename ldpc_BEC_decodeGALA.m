function c_hat = ldpc_BEC_decodeGALA(r, H, m , n, bit_node, check_node, v, u, no_reps)
%LDPC_DECODEGALA decodes the received message r passed through bsc(p)
%channel using the Gallager A message passing decoder

%m rows -> check nodes, n cols -> bit nodes of parity check matrix H
c_hat = r;
flag = 0;
for l=1:no_reps
    %1st part bit_node part
    for i=1:n %i traverse all the bit nodes
        w_i = nnz(bit_node(i,:)); %w_i is the weigth of the i-th bit node
        for j=1:w_i
            v(i,bit_node(i,j)) = c_hat(i);
        end
    end
    %2nd part check node part
    for j=1:m %%j traverse all the check nodes
        w_j = nnz(check_node(j,:)); %w_j is the weigth of the j-th check node
        sum = 0;
        erasure = zeros(1, w_j);
        for i=1:w_j
            if v(check_node(j,i), j) == 100
                erasure(i) = 1; %means that there is an erasure in the i-th edge 
            end
            
            sum = sum + v(check_node(j,i), j);
            
        end
        
        
        for i=1:w_j
            if nnz(erasure) >= 2 %means that there are more the two erasures edges in j-check node
               break;
            
            else
                if erasure(i) == 1
                    u(check_node(j,i), j) = mod(sum + v(check_node(j,i), j),2);
                end
            end
        end
    end
    
    
    for i=1:n
        
        if c_hat(i) == 100
            w_i = nnz(bit_node(i,:));
            
            for j=1:w_i
                if u(i, bit_node(i,j)) ~= 100
                   c_hat(i) = u(i, bit_node(i,j));
                       
                end
            end
        end
 
    end
    check = mod(H*c_hat.',2);
     
    if check == 0
        flag = 1;
        break;
    end
end

end
  


