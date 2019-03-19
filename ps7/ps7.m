clear
theta=.64;
alpha=1;
delta=.02;
TranMat=[.8,.2;.1,.9];
E=1;
psi=[1;0];
beta=.98;
N=2;
z=[.75;1];
%solve for equilibrium measure
mu=(eye(N)-(1-delta)*TranMat')\(E*psi);

%solve for equilibrium wage
f=@(w) sum(mu.*z.^(1/(1-theta))*(w/theta)^(theta/(theta-1)))-w;
wage=fsolve(f,1);
n=mu.*(wage./(z*theta)).^(1/(1-theta));


muPath=zeros(2,201);
nPath=zeros(2,201);
wPath=ones(1,201);
yPath=zeros(1,201);
muPath(:,1)=mu;
muPath(:,2)=[0.9*sum(mu);0]; %shock 1&2
%muPath(:,2)=[sum(mu);0]; %shock 1
%muPath(:,2)=0.9*mu; %shock 2


for j=3:201
	muPath(:,j)=E*psi+(1-delta)*TranMat'*muPath(:,j-1);
end

%solve for wage along the path
g=@(w) sum(muPath.*z.^(1/(1-theta)).*(w./theta).^(theta/(theta-1)),1)-w;
wPath=fsolve(g,wPath);

%Compute transition
nPath=(wPath./(z*theta)).^(1/(theta-1));
nTotalPath=round(sum(muPath.*nPath,1),8);
yPath=muPath.*z.*nPath.^(theta);
yTotalPath=sum(yPath,1);

%Prepare for graphs
mupath=(muPath./muPath(:,1)-1)*100;
wpath=(wPath./wPath(1)-1)*100;
mupath=mupath(:,2:201);
wpath=wpath(2:201);

t=1:200;
plot(t, mupath(1,:),t,mupath(2,:),'--')
title('Transition of measure of firms')
xlabel('Periods')
ylabel('$\%\Delta$ from initial measure of firms','Interpreter','latex')
legend({'Low Productivity','High Productivity'},'Location','southeast')
figure
plot(t,wpath)
title('Transition of wage/total production')
xlabel('Periods')
ylabel('$\%\Delta$ from initial wage/total production','Interpreter','latex')