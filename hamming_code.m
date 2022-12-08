clear;
clc;

%Parity-check and generator matrices for Hamming code (7,4)
m1 = 3;
[h1,g1,n1,k1] = hammgen(m1);
G1 = zeros(size(g1,1),size(g1,2));
H1 = zeros(size(h1,1),size(h1,2));

%make g1 in systematic form G1
for i =1:size(G1,2)
    if i > size(G1,2) - size(G1,1)
        G1(:,i - m1) = g1(:,i);
    else
        G1(:,i+m1+1) = g1(:,i);
    end
end

%make h1 in systematic form H1
for i =1:size(H1,2)
    if i > m1
        H1(:,i - m1) = h1(:,i);
    else
        H1(:,i + m1 + 1) = h1(:,i);
    end
end

%Parity-check and generator matrices for Hamming code (15,11)
m2 = 4;
[h2,g2,n2,k2] = hammgen(m2);
G2 = zeros(size(g2,1),size(g2,2));
H2 = zeros(size(h2,1),size(h2,2));

%make g2 in systematic form G2
for i =1:size(G2,2)
    if i > size(G2,2) - size(G2,1) 
        G2(:,i - m2) = g2(:,i);
    else
        G2(:,i + size(G2,1)) = g2(:,i);
    end
end

%make h2 in systematic form H2
for i =1:size(H2,2)
    if i > m2
        H2(:,i - m2) = h2(:,i);
    else
        H2(:,i + size(H2,2) - m2) = h2(:,i);
    end
end


%code Rythm
R1 = k1/n1;
R2 = k2/n2;

message_block = 1000;

c1 = zeros(message_block,n1);
c2 = zeros(message_block,n2);

for i =1:message_block
    %messages
    msg1 = zeros(1,k1);
    msg2 = zeros(1,k2);

    %encoding
    c1(i,:) = mod(msg1 * G1,2);
    c2(i,:) = mod(msg2 * G2,2);  
end

r1 = zeros(message_block,n1);
r2 = zeros(message_block,n2);

pE = linspace(0.01,1,100);
error1 = zeros(1,length(pE));
error2 = zeros(1,length(pE));

%passing messages through binary symmetric channel with error probability p
for p=1:length(pE)    
    
    
    for i = 1:message_block
        r1(i,:) = bsc(c1(i,:), pE(p));
        r2(i,:) = bsc(c2(i,:), pE(p));
        
        syn1 = mod(r1(i,:) * transpose(H1),2);
        errorbit1 = 0;
        
        syn2 = mod(r2(i,:) * transpose(H2),2);
        errorbit2 = 0;
        
        for j = 1:size(H1,2)
            if H1(:,j) == transpose(syn1)
                errorbit1 = j;
            end
        end
        for j = 1:size(H2,2)
            if H2(:,j) == transpose(syn2)
                errorbit2 = j;
            end
        end
        
        if errorbit1 ~=0
            if r1(i,errorbit1) == 0
                r1(i,errorbit1) = 1;
            else
                r1(i,errorbit1) = 0;
            end
        end
        
        if errorbit2 ~=0
            if r2(i,errorbit2) == 0
                r2(i,errorbit2) = 1;
            else
                r2(i,errorbit2) = 0;
            end
        end
        tmpr1 = r1(i,:);
        tmpc1 = c1(i,:);
        if sum(tmpr1 + tmpc1)~=0
           
            error1(p) = error1(p) + 1;
           
        end
        tmpr2 = r2(i,:);
        tmpc2 = c2(i,:);
        if sum(tmpr2 + tmpc2)~=0
            
            error2(p) = error2(p) + 1;
           
        end        
    end
    error1(p) = error1(p)/message_block;
    error2(p) = error2(p)/message_block;
end
figure(1);
plot(pE,error1);hold on;
plot(pE,error2);
set(gca,'YScale','log');
