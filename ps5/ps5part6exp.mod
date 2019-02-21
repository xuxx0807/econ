var c, k, n, z, u, y, w; 
varexo e f; 
parameters beta, alpha, theta, gam, mu, phi; 

alpha = 0.36; 
mu = 0.95; 
beta = 0.99;  
gam =0; 
theta = 2.95; 
phi = 2;

model;
beta*exp(c)^(-1)*(1-exp(u)^phi/(phi*exp(w))+alpha*exp(z)*exp(u)^alpha*exp(k)^(alpha-1)*exp(n)^(1-alpha))=exp(c(-1))^(-1);
theta*(1-exp(n))^(-gam)=exp(c)^(-1)*(1-alpha)*exp(z)*exp(u)^alpha*exp(k)^alpha*exp(n)^(-alpha);
alpha*exp(z)*exp(u)^(alpha-phi)*exp(k)^(alpha-1)*exp(n)^(1-alpha)=1;
exp(k(+1))=exp(z)*exp(u)^alpha*exp(k)^alpha*exp(n)^(1-alpha)+(1-exp(u)^phi/(phi*exp(w)))*exp(k)-exp(c);
z = mu*z(-1)+ e; 
exp(y) = exp(z)*exp(u)^alpha*exp(k)^alpha*exp(n)^(1-alpha);
w = mu*w(-1)+ f;
end; 

initval;
c=log(0.8036);
k=log(11.0836);    
n=log(0.2918);
u=log(.5);
z=0;
w=0;
e=0;
f=0;
y = log(1);
end;

shocks;
var e = 0.009^2;
var f = 0.009^2;
corr e, f = .3;
end;

steady;
stoch_simul(periods=10100); 
