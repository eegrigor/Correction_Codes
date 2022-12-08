clear;
clc;


e = 0.4;
rmin =4;
rmax = 5;
ravg_arr = rmin:0.1:rmax;
fval_max =  -Inf;
best_ravg = 0;
for z = 1:length(ravg_arr)
    ravg = ravg_arr(z);
    r = floor(ravg);
    coeff = r*(r+1-ravg)/ravg; 
    rho = @(x)coeff*x^(r-1) + (1-coeff)*x^r;
    lmax = 15;
    lmin = 2;

    %f
    f = 1./linspace(lmin, lmax, lmax-lmin+1);
    %C1:
    Aeq = ones(1, length(f));
    beq = 1;
    %C4:
    x = 0.01:0.01:e;
    %C2:
    A = zeros(length(x), length(f));
    b = ones(1,length(x));
    for i=1:length(x)
        for j=1:length(f)
            A(i,j) = e*(1-rho(1-x(i)))^j;
        end
        b(i) = x(i);
    end
    %C3:
    lb = zeros(length(f),1);
    ub = ones(length(f),1);


    [l_tmp,fval] = linprog(-f,A,b,Aeq,beq,lb,ub);
    if abs(fval) > fval_max
        fval_max = abs(fval);
        l = l_tmp;
        best_ravg = ravg;
    end
    
end

l = round(l,4); %round to fourth decimal point
num_l = zeros(1, length(l));
denom_l = 0;
for i=1:length(l)
    denom_l = denom_l + l(i)/(i+1);
end
denom_l = round(denom_l,4);

for i=1:length(l)
    num_l(i) = l(i)/((i+1));
end
num_l = round(num_l,4);

denom_l = int16(denom_l/0.0001);
num_l = int16(num_l/0.0001);

Ds = zeros(1, length(l)+2); %plus 2 because we have two rho values

for i=1:length(l)
    tmp_denom = denom_l;
    tmp_num = num_l(i);
    if tmp_num~=0
        while gcd(tmp_num, tmp_denom) ~= 1 
            div = gcd(tmp_num, tmp_denom);
            tmp_num = tmp_num/div;
            tmp_denom = tmp_denom/div;

        end
        Ds(i) = tmp_denom;
    else
        Ds(i) = 1;
    end
    
end
denom_rho = denom_l;
num_rho = zeros(1, 2);

num_rho(1) = (rmax*(rmax+1-best_ravg)/best_ravg)/rmin;
num_rho(2) = (1-rmax*(rmax+1-best_ravg)/best_ravg)/rmax;
num_rho = round(num_rho,4);
num_rho = int16(num_rho/0.0001);
for i=1:length(num_rho)
    tmp_denom = denom_rho;
    tmp_num = num_rho(i);
    if tmp_num~=0
        while gcd(tmp_num, tmp_denom) ~= 1 
            div = gcd(tmp_num, tmp_denom);
            tmp_num = tmp_num/div;
            tmp_denom = tmp_denom/div;

        end
        Ds(end-mod(i,2)) = tmp_denom;
    else
        Ds(end-mod(i,2)) = 1;
    end
    
end

n = lcms(Ds);
