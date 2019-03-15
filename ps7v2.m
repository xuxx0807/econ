%1.1
E = 1;
psi = [1, 0];
delta = .02;
F = [.8, .2; .1, .9];
I = eye(2);

mustar = E*psi*inv(I-(1-delta)*F)

%1.2
Z=[.75;1];
alpha = 1;
theta = .64;
nstar = zeros(2,1);
ystar = zeros(2,1);
pistar = zeros(2,1);

wold = 1E-15;
wnew = 1;
step = 1;
wd = 0;

while step > 1E-15
    for i =1:2
        nstar(i) = (wnew/(Z(i)*theta))^(1/(theta-1));
        ystar(i) = Z(i)*nstar(i)^theta;
        pistar(i) = ystar(i)- wnew*nstar(i);
    end
    wd = alpha*sum(mustar*ystar);
    if wd<wnew
        wnew=wold;
        step = step/2;
    else
        wold = wnew;
        wnew = wnew+step;
    end
end

wnew

%2
%mushock = [.9*sum(mustar), 0];
%3
%mushock = [sum(mustar), 0];
mushock = .9*mustar;
mupath = zeros(1,2,201);
mupath (:,:,1) = mushock;

for i=1:200
    mupath(:,:,i+1) = E*psi+ mupath(:,:,i)*(1-delta)*F;
end

wpath = zeros(1,201);
yspath = zeros(2,1,201);
nspath = zeros(2,1,201);
ypath = zeros(1,201);
npath = zeros(1,201);

for i=1:201
    nstar2 = zeros(2,1);
    ystar2 = zeros(2,1);
    pistar2 = zeros(2,1);
    wold2 = 1E-15;
    wnew2 = 1;
    step2 = 1;
    wd2 = 0;
    while step2 > 1E-15
        for j =1:2
            nstar2(j) = (wnew2/(Z(j)*theta))^(1/(theta-1));
            ystar2(j) = Z(j)*nstar2(j)^theta;
            pistar2(j) = ystar2(j)- wnew*nstar2(j);
        end
        wd2 = alpha*sum(mupath(:,:,i)*ystar2);
        if wd2<wnew2
            wnew2=wold2;
            step2 = step2/2;
        else
            wold2 = wnew2;
            wnew2 = wnew2+step2;
        end
    end
    wpath(1,i)=wnew2;
    yspath(:,:,i) = ystar2;
    nspath(:,:,i) = nstar2;
    ypath(1,i) = sum(mupath(:,:,i)*yspath(:,:,i));
    npath(1,i) = sum(mupath(:,:,i)*nspath(:,:,i));
end
x = zeros(1,201);
for i=1:201;
    x(i)=i;
end

normwpath = wpath/wnew;
normnpath = npath/sum(mustar*nstar);
normypath = ypath/sum(mustar*ystar);
%plot(x,wpath,'g',x,npath,'b',x,ypath,'r')

%ypathboth = ypath;
%ypath1 = ypath;
%ypath2 = ypath;

%normypathboth = normypath;
%normypath1 = normypath;
%normypath2 = normypath;



plot(x,ypathboth,'g',x,ypath1,'b',x,ypath2,'y')

