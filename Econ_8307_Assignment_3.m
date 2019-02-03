%Q1
N=10;
epsilon=1;
delta=.1;
beta=1;
gamma=.1;
z=1:N;
transitMatrix=eye(N);
alpha=zeros(N,1);
phi=ones(N,1)/N;
mu=zeros(N,1);

for i=1:N
	for j=1:N
		alpha(i)=alpha(i)+max(0,beta-gamma*(z(j)-z(i))^2);
	end
end
for i=1:N
	for j=1:N
		transitMatrix(i,j)=max(0,(beta-gamma*(z(j)-z(i))^2)/alpha(i))
	end
end

mu=((1-delta)*transitMatrix.'-eye(N))\(-epsilon*phi)

%Q3
clear
N=5;
w=1;
alpha=.7;
beta=.95;
E=1;
lambda=.1;
tau=0;
z=1:N
phi=ones(N,1)/N;
transitMatrix=0.05*ones(N)+0.75*eye(N);
muHat=((1-lambda)*transitMatrix.'-eye(N))\(-phi)

n=@(z) (1./(alpha*z)).^(1/(alpha-1));
labor=n(z)*phi;
