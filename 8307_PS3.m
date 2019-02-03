%Macro III
%Assignment 3
%Fred Xu and Jonah Coste
%2/1/19

%Question 1
%Part 4
N = 10;
e = 1;
delta = .1;
beta = 1;
gamma = .1;

psi = zeros(1,N);
T = zeros(N);
I = eye(N);
A= zeros(1,N);

relprob=@(new, old) max(0, beta - gamma*(new-old)^2);

for i = 1:N
    psi(i) = 1/N;
    for j= 1:N
        A(i) = A(i) + relprob(j, i);
    end
end

for i=1:N
    for j = 1:N
        T(i,j) = relprob(j, i) / A(i);
    end
end

mustar = e*psi*inv(I-(1-delta)*T);

mustar'

%Question 3
%Part 1
N = 5;
p = .8;
w = 1;
alpha = .7;
beta = .95;
lambda = .1;

Z = zeros(1,N);
nstar = zeros(N,1);
for i =1:N
    Z(i) = i;
    nstar(i) = (w/(Z(i)*alpha))^(1/(alpha-1));
end

nstar

T = zeros(N);
T=T+(1-p)/(N-1)+eye(N)*(p-(1-p)/(N-1));
I=eye(N);
pi = zeros(N,1);
for i =1:N
     pi(i)=Z(i)*nstar(i)^alpha - w*nstar(i);
 end
V = inv(I-beta*(1-lambda)*T)*pi
