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

%Q3
clear
N=5;
w=1;
alpha=.7;
beta=.95;
E=1;
lambda=.1;
tau=0.5;
z=1:N;
phi=ones(N,1)/N;
transitMatrix=0.05*ones(N)+0.75*eye(N);
muHat=((1-lambda)*transitMatrix.'-eye(N))\(-phi);

n=@(z) (1./(alpha*z)).^(1/(alpha-1));
labor=n(z)*phi;

indexat=@(fun,index) fun(:,index);
freeEntry=@(wage) beta*indexat(valueIteration(wage),1).'*phi-E;
wStar=fsolve(freeEntry,0)
value=valueIteration(wStar);
nGridNum=100;
nGrid=0:.5/(nGridNum-1):.5;

mu=zeros(N,nGridNum);
muNew=zeros(N,nGridNum);
epsilon=100;
iteration=0;
while epsilon>1e-5&&iteration<500
	for a=1:N
		for b=1:nGridNum
			for i=1:N
				for j=1:nGridNum
					muNew(a,b)=muNew(a,b)+(1-lambda)*(value(a,j+100)==nGrid(b))*transitMatrix(i,a)*mu(i,j);
				end
			end
			muNew(a,b)=muNew(a,b)+phi(a)*(b==1);
		end
	end
	epsilon=sum(sum(abs(muNew-mu)));
	mu=muNew;
	iteration=iteration+1;
end

nTotal=sum(sum(mu.*value(:,101:200)));
nLoss=mu.*value(:,101:200)*lambda+mu.*value(:,101:200)*(1-lambda).*max(0,nGrid-value(:,101:200));


value1=valueIteration(1);
function valMat = valueIteration (w)
N=5;
alpha=.7;
beta=.95;
lambda=.1;
tau=0;
z=1:N;
phi=ones(N,1)/N;
transitMatrix=0.05*ones(N)+0.75*eye(N);
nGridNum=100;
nGrid=0:.5/(nGridNum-1):.5;
epsilon=100;
valueMatrix=zeros(N,nGridNum);
iteration=0;
ind=zeros(N,nGridNum);
while (epsilon>1e-5) && (iteration<500)
for i=1:N
	for j=1:nGridNum
		[newValueMatrix(i,j),ind(i,j)]=max(z(i)*nGrid.^alpha-w*nGrid-tau*w*max([zeros(1,nGridNum);nGrid(j)-nGrid],[],1)+beta*((1-lambda)*transitMatrix(i,:)*valueMatrix-lambda*tau*w*nGrid));
	end
end
epsilon=sum(sum(abs(newValueMatrix-valueMatrix)));
valueMatrix=newValueMatrix;
iteration=iteration+1;
end
nMatrix=nGrid(ind);
valMat=[valueMatrix,nMatrix];
end

