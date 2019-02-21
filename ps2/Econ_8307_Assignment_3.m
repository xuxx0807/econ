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
		transitMatrix(i,j)=max(0,(beta-gamma*(z(j)-z(i))^2)/alpha(i));
	end
end

mu=((1-delta)*transitMatrix.'-eye(N))\(-epsilon*phi);
muStar=epsilon/delta.*mu;

%Q3
clearvars -except muStar;
global N alpha beta E transitMatrix omega z lambda phi nGridNum nmin nmax;
nGridNum=100;
nmin=0;
nmax=0.5;
omega=1;
N=5;
w=1;
alpha=.7;
beta=.95;
E=1;
lambda=.1;
z=1:N;
phi=ones(N,1)/N;
transitMatrix=0.05*ones(N)+0.75*eye(N);
muHat=((1-lambda)*transitMatrix.'-eye(N))\(-phi);

n=@(z) (1./(alpha*z)).^(1/(alpha-1));
labor=n(z)*phi;


eqm1=freeEntry(0);
%eqm 1st element is wage, 2nd is value matrix, 3rd is labor matrix
mass1=measure(eqm1);
desRate1=job(eqm1);
employ1=hh(eqm1);
util1=log(eqm1{1}/omega)-employ1;

eqm2=freeEntry(0.5);
mass2=measure(eqm2);
desRate2=job(eqm2);
employ2=hh(eqm2);
util2=log(eqm2{1}/omega)-employ2;
percentage2=exp(omega*(employ2-employ1))*eqm1{1}/eqm2{1};

eqm3=freeEntry(1);
mass3=measure(eqm3);
desRate3=job(eqm3);
employ3=hh(eqm3);
util3=log(eqm3{1}/omega)-employ3;
percentage3=exp(omega*(employ3-employ1))*eqm1{1}/eqm3{1};


function [valueMatrix, nMatrix] = valueIteration (w,t)
%Given wage and firing cost rate, return matrices of value and optima; n_t given z_t and n_{t-1}, using value function iteration
global N alpha beta lambda z phi transitMatrix nGridNum nmin nmax;
nGrid=nmin:nmax/(nGridNum-1):nmax;
epsilon=100;
valueMatrix=zeros(N,nGridNum);
iteration=0;
ind=zeros(N,nGridNum);
while (epsilon>1e-5) && (iteration<500)
for i=1:N
	for j=1:nGridNum
		[newValueMatrix(i,j),ind(i,j)]=max(z(i)*nGrid.^alpha-w*nGrid-t*w*max([zeros(1,nGridNum);nGrid(j)-nGrid],[],1)+beta*((1-lambda)*transitMatrix(i,:)*valueMatrix-lambda*t*w*nGrid));
	end
end
epsilon=sum(sum(abs(newValueMatrix-valueMatrix)));
valueMatrix=newValueMatrix;
iteration=iteration+1;
end
nMatrix=nGrid(ind);
end

function eqm = freeEntry (t)
%Given firing cost rate, return equilibrium wage, value matrix, and labor matrix that satisfying free entry
global N phi E beta;
indexat=@(fun,index) fun(:,index);
f=@(wage) beta*indexat(valueIteration(wage,t),1).'*phi-E;
wStar=fsolve(f,0);
[valueMatrix, nMatrix]=valueIteration(wStar,t);
eqm={wStar,valueMatrix,nMatrix};
end

function mu = measure (eqm)
%Given equilibrium satisfying free entry, return distribution over firm types, using iteration
global N nGridNum lambda phi transitMatrix nmin nmax;
nGrid=nmin:nmax/(nGridNum-1):nmax;
mu=zeros(N,nGridNum);
muNew=zeros(N,nGridNum);
epsilon=100;
iteration=0;
value=eqm{3};
while epsilon>1e-5&&iteration<500
	for a=1:N
		for b=1:nGridNum
			for i=1:N
				for j=1:nGridNum
					muNew(a,b)=muNew(a,b)+(1-lambda)*(value(i,j)==nGrid(b))*transitMatrix(i,a)*mu(i,j);
				end
			end
			muNew(a,b)=muNew(a,b)+phi(a)*(b==1);
		end
	end
	epsilon=sum(sum(abs(muNew-mu)));
	mu=muNew;
	muNew=zeros(N,nGridNum);
	iteration=iteration+1;
end
end

function desRate = job (eqm)
%Given equilibrium satisfying free entry, return job destruction rate
global lambda nGridNum nmin nmax N transitMatrix;
nGrid=nmin:nmax/(nGridNum-1):nmax;
mass=measure(eqm);
value=eqm{3};
nTotal=sum(sum(mass.*value));
nmat=repmat(nGrid,N,1);
%nLoss=sum(sum(mass.*value*lambda+mass.*max(0,nGrid-value)));
nLoss=sum(sum(mass.*eqm{3}*lambda+mass.*max(0,nmat-value)));
desRate=nLoss/nTotal;
end

function employ = hh (eqm)
%Given equilibrium satisfying free entry, return employment, using labor market clearing condition
global omega E alpha N z nGridNum
mu=measure(eqm);
value=eqm{3};
mat=repmat(z,nGridNum,1).';
laborClear=@(M) eqm{1}/omega+M*E-M*sum(sum(mat.*(value.^alpha).*mu));
mStar=fsolve(laborClear,0);
employ=sum(sum(mStar*mu.*eqm{3}));
end