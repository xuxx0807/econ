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
beta*c^(-1)*(1-u^phi/(phi*exp(w))+alpha*exp(z)*u^alpha*k^(alpha-1)*n^(1-alpha))=c(-1)^(-1);
theta*(1-n)^(-gam)=c^(-1)*(1-alpha)*exp(z)*u^alpha*k^alpha*n^(-alpha);
alpha*exp(z)*u^(alpha-phi)*k^(alpha-1)*n^(1-alpha)=1;
k(+1)=exp(z)*u^alpha*k^alpha*n^(1-alpha)+(1-u^phi/(phi*exp(w)))*k-c;
z = mu*z(-1)+ e; 
y = exp(z)*u^alpha*k^alpha*n^(1-alpha);
w = mu*w(-1)+ f;
end; 

initval;
c=0.8036;
k=11.0836;    
n=0.2918;
u=.5;
z=0;
w=0;
e=0;
f=0;
y = 1;
end;

shocks;
var e = 0.009^2;
var f = 0.009^2;
corr e, f = .3;
end;

steady;
stoch_simul(periods=10100); 
