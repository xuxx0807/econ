%Macro III
%Assignment 2
%Fred Xu and Jonah Coste
%1/23/19

%Question 1
beta = .97;
delta = .1;
theta = .3;
foc=@(policy) [policy(2)^theta-delta*policy(2)-policy(1);
	policy(1)*(beta*(theta * policy(2)^(theta-1))-1)+1-delta];
%First order conditions in the steady state, taking vector policy as argument
%consumption and capital are the first and second element of the policy vector
initial_guess=[1;1];
steady_state=fsolve(foc,initial_guess)

%Question 2
%2.1
%V^e = u(w_h) + \beta [\pi_e V^e + (1-\pi_e) V^u]
%V^u = u(w_l) + \beta [\pi_u V^u + (1-\pi_u) V^e]

%2.2
u=@(c) log(c);
beta=.96;
wh=1;
wl=.9;
prob_e=.95;
prob_u=.9;
coef=[beta*prob_e-1, beta*(1-prob_e); beta*(1-prob_u), beta*prob_u-1];
%Left hand side of the linear equations whose first variable is value when employed
%and second is value when unemployed
const=[-u(wh);-u(wl)];
v=inv(coef)*const

%2.3
prob_matrix=[prob_e, 1-prob_u; 1-prob_e, prob_u];
[V,D]=eig(prob_matrix);
unemploymentRate=V(2,1)/(V(1,1)+V(2,1))

%Question 3
%n=@(z) (w./(alpha*z)).^(1/(alpha-1));
%value=z.*n(z).^alpha-w*n(z)+beta*(1-lambda)*probMat*value;

%Question 4
%4.1
N=5;
z=zeros(N,1);
for i=1:N
	z(i)=i/N;
end
probMat=zeros(N)+.05+eye(N)*(.8-.05);
w=1;
alpha=.7;
beta=.95;
lambda=.1;
coef=beta*(1-lambda)*probMat-eye(N);
n=@(z) (w./(alpha*z)).^(1/(alpha-1));
const=-(z.*n(z).^alpha-w*n(z));
vStar=inv(coef)*const

%4.2
t=1;
T=100;
value=zeros(N,T);
L=zeros(T,1);
D=zeros(T-1,1);
indL=0;
indD=0;
epsilon=1e-4;
while t<=T
	value(:,t+1)=z.*n(z).^alpha-w*n(z)+beta*(1-lambda)*probMat*value(:,t);
	L(t)=sum(abs(value(:,t)-vStar));
	D(t)=sum(abs(value(:,t+1)-value(:,t)));
	if D(t)<epsilon
		indD=indD+1;
	end
	if L(t)<epsilon
		indL=indL+1;
	end
	if indL==1
		L_iter=t
	end
	if indD==1
		D_iter=t
	end
	t=t+1;
end
plot(1:T,L,1:T-1,D(1:T-1))