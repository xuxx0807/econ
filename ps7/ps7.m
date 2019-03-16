global theta alpha delta psi beta E;
theta=.64;
alpha=1;
delta=.02;
TranMat=[.8,.2;.1,.9];
E=1;
psi=[1;0];
beta=.98;
N=2;
z=[.75;1];
mu=(eye(N)-(1-delta)*TranMat')\(E*psi);
laborMarket = @(w) sum(mu.*(w/(theta.*z).^(1/(theta-1))),'all')-alpha/theta;
guess=1;
wage=fsolve(laborMarket,guess);

muPath=zeros(2,200);
nPath=zeros(2,200);
wPath=zeros(1,200);
yPath=zeros(1,200);
muPath(:,1)=[45;0];
for j=2:200
	muPath(:,j)=E*psi+(1-delta)*TranMat'*muPath(:,j-1);
end