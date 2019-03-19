%Parameters
gamma = 1;
alpha = .36;
theta = .5;
beta = .95;
Z = [.5, 1];
T = [.9, .1; .1, .9];
wguess = 1;
omega = 1;


%Grid settings
gs = 10;
gmax = 2;

gv = zeros(1,gs+1);
for i = 1:gs+1
    gv(1,i)= gmax*(i-1)/gs;  
end

phi = zeros(2, gs+1, gs+1, gs+1);
phi(:,1,1,1) = [.5;.5];
%Find equilibrium wage using free entry condition
w = wguess*2;
wcons = 1E-15;
%wcons = wguess;
%w = wguess+1;

Vguess = zeros(2, gs+1, gs+1, gs+1);
Nstari = zeros(2, gs+1, gs+1, gs+1);
Vnew = zeros(2, gs+1, gs+1, gs+1);
muold = zeros(2, gs+1, gs+1, gs+1);
muold(1, 1, 1, 1) = .5;
muold(2, 1, 1, 1) = .5;
munew = muold;

%iterates to find wage that satisfies ls = ld
while abs(log(w/wcons)) > .02
w = (w+wcons)/2;
%w = wcons

loss = 100;

%find value function given wage 
while loss > 1
        for ti=1:gs+1
            for ni=1:gs+1
                for hi=1:gs+1
                    temp = zeros(2, gs+1);
                    for newni=1:gs+1
                        newn = gv(newni);
                        newt = min(newn, gv(ti) + theta*gv(hi));
                        newh = max(newn - gv(ni), 0);
                        [dt, newti] = min(abs(gv-newt));
                        [dh, newhi] = min(abs(gv-newh));
                        temp(:, newni)= Z.'*(newn+gamma*newt)^alpha - w*newn + beta*T*Vguess(:,newti,newni,newhi);
                    end
                    [maxV, maxVni] = max(temp.');
                    Nstari(:, ti, ni, hi) = maxVni;
                    Vnew(:, ti, ni, hi) = maxV;
                end
            end
        end
    loss = sum(abs(Vnew-Vguess), 'all');
    Vguess = Vnew;
end

%find ss mu given value function
loss2 = 1;
while loss2 > .00001
muold = munew;
munew = zeros(2, gs+1, gs+1, gs+1);
for toldi=1:gs+1
    for noldi=1:gs+1
        for holdi=1:gs+1
            for zoldi=1:2
                told = gv(toldi);
                nold = gv(noldi);
                hold = gv(holdi);
                nnewi = Nstari(zoldi, toldi, noldi, holdi);
                nnew = gv(nnewi);
                hnew = max(nnew-nold,0);
                [dh, hnewi] = min(abs(gv-hnew));
                tnew = min(nnew, told+theta*hold);
                [dt, tnewi] = min(abs(gv-tnew));
                munew(:, tnewi, nnewi, hnewi) = munew(:, tnewi, nnewi, hnewi) + T(zoldi,:).'*muold(zoldi, toldi, noldi, holdi);
            end
        end
    end
end
loss2 = sum(abs(munew-muold), 'all');
end

%find implied W from consumer problem
Ystar = zeros(2, gs+1, gs+1, gs+1);
for toi=1:gs+1
    for noi=1:gs+1
        for hoi=1:gs+1
            for zi=1:2
                to = gv(toi);
                ho = gv(hoi);
                nn = gv(Nstari(zi, toi, noi, hoi));
                tn = min(nn, to+theta*ho);
                Ystar(zi, toi, noi, hoi) = Z(zi)*(nn+gamma*tn)^alpha;
            end
        end
    end
end
Y = sum(munew.*Ystar, 'all');
wcons = omega*Y;
end

maxchosen = max(Nstari, [], 'all'); %check! if =11 then increase gmax
V= Vnew;
mustar = munew;
Emp = sum(mustar.*gv(Nstari), 'all');

