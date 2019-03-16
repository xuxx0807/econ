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
f=@(w) sum(mu.*z.^(1/(1-theta))*(w/theta)^(theta/(theta-1)))-w;
wage=fsolve(f,1);
for j=2:200
	muTran(:,j)=E*psi+(1-delta)*TranMat'*muTran(:,j-1);
end