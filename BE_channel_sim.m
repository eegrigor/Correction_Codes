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

%spy(H);
% decode the message r using Gallager A message passing decoder
check_node = []; %stores the connected bit nodes to the j-th check node
n = n*z;
m = m*z;
check_weights = zeros(1,m);

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
v = 100*ones(n,m);
u = 100*ones(n,m);
no_reps = 10; %number of iterations
%%%%%%%%%

%code rate
R = (n-m)/n;

block_length = 100;
eps = linspace(0.01, 0.99, 99);
no_reps = [10 30 100];
ber = zeros(length(no_reps), length(eps));
for l=1:length(no_reps)
    for pE=1:length(eps)
        errors = 0;
        for i=1:block_length
            noise = randsrc(1, n, [1 0; eps(pE) 1-eps(pE)]);
            % we send the all zero codeword
            c = zeros(1,n);
            %pass the codeword c through a BEC(eps) channel with eps probability of
            %erasure
            %received
            r = c;
            for j=1:length(noise)
                if noise(j) == 1
                    r(j) = 100; % this models the erasure
                end
            end 

            c_hat = ldpc_BEC_decodeGALA(r, H, m , n, bit_node, check_node, v, u, no_reps(l));

            errors = errors + nnz(c_hat);
        end
        ber(l,pE) = errors/(block_length*n);
    end
end
figure(1);
plot(eps, ber(1,:));hold on;
set(gca, 'YScale', 'log')
plot(eps, ber(2,:));hold on;
set(gca, 'YScale', 'log')
plot(eps, ber(3,:));
set(gca, 'YScale', 'log')
xlabel('Erasure Probability');
ylabel('BER');
legend('reps = 10', 'reps = 30', 'reps = 100');

