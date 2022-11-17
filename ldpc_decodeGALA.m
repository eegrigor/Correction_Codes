function c_hat = ldpc_decodeGALA(r, h, m , n)
%LDPC_DECODEGALA decodes the received message r passed through bsc(p)
%channel using the Gallager A message passing decoder

%m rows -> check nodes, n cols -> bit nodes of parity check matrix H

%initialize 2d check node matrix
check_node = [];
check_weights = zeros(1,n);

for i=1:size(h,1)
    temp = h(i,1);
    check_weights(temp) = check_weights(temp)+1;
    check_node(temp, check_weights(temp)) = h(i,2);
end

%initialize 2d bit node matrix
bit_node = [];
bit_weights = zeros(1,n);

for i=1:size(h,1)
    temp = h(i,2);
    bit_weights(temp) = bit_weights(temp)+1;
    bit_node(temp, bit_weights(temp)) = h(i,1);
end

%for every bit node initialize matrices v and u with length = #check nodes
v = zeros(n,m);
u = zeros(n,m);
l = 1; %number of iterations
for l=1:l
    %1st part bit_node part
    for i=1:n
        w_i = nnz(bit_node(i,:));
        if l==1
            for j=1:w_i
                v(i,bit_node(i,j)) = r(i);
            end
            
        else
            %check if u(i) are equal
            temp = u(i,bit_node(i,1));
            is_equal = 1; %0 -> no, 1 -> yes
            for j=2:w_i
                if u(i,bit_node(i,j)) ~= temp
                    is_equal = 0;
                end
                    
            end
            for j=1:w_i
                if is_equal == 1
                    v(i,bit_node(i,j)) = u(i,bit_node(i,j));
                else 
                    v(i,bit_node(i,j)) = r(i);
                end
            end
            
        end
    end
    %2nd part check node part
    for i=1:m
        w_j = nnz(check_node(i,:));
        sum = 0;
        for j=1:w_j
            sum = sum + v(check_node(i,j), i);
        end
        for j=1:w_j
            u(check_node(i,j), i) = mod(sum + v(check_node(i,j), i),2);
        end
    end
end
c_hat = zeros(1,n);

for i=1:n
    w_i = nnz(bit_node(i,:));
    no_ones = 0;
    no_zeros = 0;
    for j=1:w_i
        if u(i, bit_node(i,j)) == 0
           no_zeros = no_zeros + 1;
        elseif u(i, bit_node(i,j)) == 1
           no_ones = no_ones + 1;
        end
        
    end
    if no_zeros > no_ones
        c_hat(i) = 0;
    elseif no_zeros < no_ones
        c_hat(i) = 1;
    else
        c_hat(i) = r(i);
    end
        
end


end
  


