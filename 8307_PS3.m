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

%Part 2 a
T = zeros(N);
T=T+(1-p)/(N-1)+eye(N)*(p-(1-p)/(N-1));
I=eye(N);
psi= zeros(1,N);
for i =1:N
    psi(i) = 1/N;
end

mustar = psi*inv(I-(1-lambda)*T);
mustar1 = mustar/sum(mustar);
mustar1'

%Part 2 b
%First find endogenous wage with entry condition.
E = 1;
wguess = 1;
ECalc = 0;
step = 1;
pi = zeros(N,1);

while step > 1E-15
    for i =1:N
        nstar(i) = (wguess/(Z(i)*alpha))^(1/(alpha-1));
        pi(i)=Z(i)*nstar(i)^alpha - wguess*nstar(i);
    end
    V = inv(I-beta*(1-lambda)*T)*pi;
    Ecalc = beta*psi*V;
    if Ecalc-E >0
        wguess = wguess + step;
    else
        wguess = wguess - step;
        step = step/2;
    end
end
w = wguess;

%Now find steady state distribution.

%Find optimal n matrix
gridsize = 200;
gridmax = .2;
nvalues = zeros(gridsize+1,1);
for i = 1:gridsize+1
    nvalues(i)= gridmax*(i-1)/gridsize; 
end
nstar2 = zeros(N, gridsize+1);
temp = zeros(gridsize+1,1);
for i = 1:N
    for j = 1:gridsize+1
        for k = 1:gridsize+1
            temp(k) = Z(i)*nvalues(k)^alpha - w*nvalues(k);
        end
        [M, I]= max(temp);
        nstar2(i,j) = nvalues(I);
    end
end

%iterate to find steady state measure
mmold = ones(N, gridsize+1);
loss=1;

while loss > 1E-15
mmnew = zeros(N, gridsize+1);
for i = 1:N
    for j = 1:gridsize+1
        for ii = 1:N
            for jj = 1:gridsize+1
                if nstar2(ii,jj) == nvalues(j)
                    mmnew(i,j) = mmnew(i,j) + (1-lambda)*T(i,ii)*mmold(ii,jj);
                end
            end
        end
        if nvalues(j) == 0
            mmnew(i,j) = mmnew(i,j) + psi(i);
        end
    end
end
loss = sum(abs(mmnew-mmold), 'all');
mmold = mmnew;
end

mm = mmnew/sum(mmnew,'all');
mu = sum(mm')

jobs = 0;
dest = 0;
for i= 1:N
    for j = 1:gridsize+1
        jobs = jobs + mm(i,j)*nstar(i);
        for k = 1:N
            if k < i
                dest = dest + mm(i,j)*T(i,k)*(nstar(i)-nstar(k));
            end
        end
    end
end

propdest = dest/jobs


