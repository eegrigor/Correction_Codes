clear;
clc;
%Parity check matrix in the base matrix form
load base_matrices/NR_1_2_20.txt;
B = NR_1_2_20;
z = 20; % expansion factor
H = base_to_parity_check(B,z);
[m,n] = size(B);

h = sparse(H);
[h1,h2] = find(h);
h = cat(2, h1,h2);

spy(H)


%code rate
R = (n-m)/n;

%message of len=n-m
msg = randi([0 1], 1, (n-m)*z);

%encode the message -> c using the Base Matrix B as in the 5G ldpc
%standars
c = nrldpc_encode(B,z, msg); 

%pass the codeword c through a BSC(p) channel with p probability of flip
p=0.04;

r = bsc(c, p);
cr = mod(c+r,2);
nnz(cr)
%%%%%%%%
% decode the message r using Gallager A message passing decoder
check_node = []; %stores the connected bit nodes to the j-th check node
n = n*z;
m = m*z;
check_weights = zeros(1,n);

for i=1:size(h,1)
    temp = h(i,1);
    check_weights(temp) = check_weights(temp)+1;
    check_node(temp, check_weights(temp)) = h(i,2);
end

%initialize 2d bit node matrix
bit_node = []; %stores the connected check nodes to the i-th bit node
bit_weights = zeros(1,n);

for i=1:size(h,1)
    temp = h(i,2);
    bit_weights(temp) = bit_weights(temp)+1;
    bit_node(temp, bit_weights(temp)) = h(i,1);
end

%for every bit node initialize matrices v and u with length = #check nodes
v = 2.2*ones(n,m);
u = 2.2*ones(n,m);
no_reps = 100; %number of iterations
for l=1:no_reps
    %1st part bit_node part
    for i=1:n %i traverse all the bit nodes
        w_i = nnz(bit_node(i,:)); %w_i is the weigth of the i-th bit node
        if l==1
            for j=1:w_i
                v(i,bit_node(i,j)) = r(i);
            end
            
        else
            %check if u(i) are equal
            
            is_equal = ones(1,w_i); %0 -> no, 1 -> yes
            temp = 0;
            for j=1:w_i
                if j == w_i
                    temp = u(i,bit_node(i,1));
                else
                    temp = u(i,bit_node(i,j+1));
                end
                for t=1:w_i
                    if t == j
                        continue;
                    end
                    if u(i,bit_node(i,t)) ~= temp
                        is_equal(j) = 0;
                        break;
                    end
                end    
            end
            for j=1:w_i
                if is_equal(j) == 1
                    if j == w_i
                        v(i,bit_node(i,j)) = u(i,bit_node(i,1));
                    elseif j == 1
                        v(i,bit_node(i,j)) = u(i,bit_node(i,2));
                    else
                        v(i,bit_node(i,j)) = u(i,bit_node(i,j+1));
                        
                    end
                
                else 
                    v(i,bit_node(i,j)) = r(i);
                end
            end
            
        end
    end
    %2nd part check node part
    for j=1:m %%j traverse all the check nodes
        w_j = nnz(check_node(j,:)); %w_j is the weigth of the j-th check node
        sum = 0;
        for i=1:w_j
            sum = sum + v(check_node(j,i), j);
        end
        for i=1:w_j
            u(check_node(j,i), j) = mod(sum + v(check_node(j,i), j),2);
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
% c_hat = ldpc_decodeGALA(r, h, m*z, n*z);
data_len = n - m
diff = mod(c_hat(1:data_len) + c(1:data_len),2);
nnz(diff)
