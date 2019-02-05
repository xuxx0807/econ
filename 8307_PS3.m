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

%Part 2 b and 3,4,5
N = 5;
p = .8;
alpha = .7;
beta = .95;
lambda = .1;
E = 1;
wguess = 1;
tau = 0.0;

Z = zeros(1,N);

for i =1:N
    Z(i) = i;
end

T = zeros(N);
T=T+(1-p)/(N-1)+eye(N)*(p-(1-p)/(N-1));
I=eye(N);
psi= zeros(1,N);
for i =1:N
    psi(i) = 1/N;
end

%Find equilibrium wage
step = 1;
gridsize = 100;
gridmax = 100;

while step > .001
% Find Value function given wguess

nvalues = zeros(gridsize+1,1);
for i = 1:gridsize+1
    nvalues(i)= gridmax*(i-1)/gridsize; 
end

nstar = zeros(N, gridsize+1);
Vold = zeros(N, gridsize+1);
temp = zeros(gridsize+1,1);
loss = 1;

while loss > .01
Vnew = zeros(N, gridsize+1);
for i = 1:N
    for j = 1:gridsize+1
        temp = zeros(gridsize+1,1);
        for k = 1:gridsize+1
            temp(k) = Z(i)*nvalues(k)^alpha - wguess*nvalues(k) - tau*wguess*max(0,nvalues(j)-nvalues(k)) + beta*(1-lambda)*T(i,:)*Vold(:,k) - beta*lambda*tau*wguess*nvalues(k);
        end
        [M, I]= max(temp);
        nstar(i,j) = nvalues(I);
        Vnew(i,j) = Z(i)*nstar(i,j)^alpha - wguess*nstar(i,j) - tau*wguess*max(0,nvalues(j)-nstar(i,j)) + beta*(1-lambda)*T(i,:)*Vold(:,I) - beta*lambda*tau*wguess*nstar(i,j);
    end
end
loss = sum(abs(Vnew-Vold), 'all');
Vold = Vnew;
end
V = Vnew;
gridmax = round(max(nstar,[],'all')*1.1,2,'significant');
Ecalc = beta*psi*V(:,1);
if Ecalc-E >0
    wguess = wguess + step;
else
    wguess = wguess - step;
    step = step/2;
end
end

%iterate to find steady state measure
muold = ones(N, gridsize+1);
loss2=1;

while loss2 > .01
munew = zeros(N, gridsize+1);
for i = 1:N
    for j = 1:gridsize+1
        for ii = 1:N
            for jj = 1:gridsize+1
                if nstar(ii,jj) == nvalues(j)
                    munew(i,j) = munew(i,j) + (1-lambda)*T(i,ii)*muold(ii,jj);
                end
            end
        end
        if nvalues(j) == 0
            munew(i,j) = munew(i,j) + psi(i);
        end
    end
end
loss2 = sum(abs(munew-muold), 'all');
muold = munew;
end

mm = munew/sum(munew,'all');
mu = sum(mm')

jobs = 0;
dest = 0;
for i= 1:N
    for j = 1:gridsize+1
        jobs = jobs + mm(i,j)*nstar(i,j);
        dest = dest + lambda*mm(i,j)*nstar(i,j);
        for k = 1:N
                kI = find(nvalues==nstar(i,j));
                dest = dest + (1-lambda)*mm(i,j)*T(i,k) * max(0,(nstar(i,j)-nstar(k,kI)));
            end
        end
end
propdest = dest/jobs

omega = 1;
c = wguess/omega;
y = 0;
emp = 0;
for i = 1:N
    for j = 1:gridsize+1
        y = y + Z(i)*nstar(i,j)^alpha * munew(i,j);
        emp = emp + nstar(i,j)*munew(i,j);
    end
end
Mstar = c/(y-E);
emp = emp*Mstar

increase = exp((.9017+emp) - log(wguess)) -1







