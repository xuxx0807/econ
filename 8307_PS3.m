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

