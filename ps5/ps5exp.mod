% This is the file from the notes
% Remember to save files like this where MATLAB can find them

% Tell dynare what variables and parameters there are
var c, k, n, z, u, y; 
varexo e; 
parameters beta, alpha, theta, gam, mu, phi; 
 
% Give values to the parameters
alpha = 0.36; 
mu = 0.95; 
beta = 0.99;  
gam =0; 
theta = 2.95; 
phi = 4;

% The dynamic equations go here
% Note the conventions for time. State variables like capital go from t to t+1
% Endogenous variables like c and z go from t-1 to t. Labor n does not appear dynamically

model;
beta*exp(c)^(-1)*(1-exp(u)^phi/phi+alpha*exp(z)*exp(u)^alpha*exp(k)^(alpha-1)*exp(n)^(1-alpha))=exp(c(-1))^(-1);
theta*(1-exp(n))^(-gam)=exp(c)^(-1)*(1-alpha)*exp(z)*exp(u)^alpha*exp(k)^alpha*exp(n)^(-alpha);
alpha*exp(z)*exp(u)^(alpha-phi)*exp(k)^(alpha-1)*exp(n)^(1-alpha)=1;
exp(k(+1))=exp(z)*exp(u)^alpha*exp(k)^alpha*exp(n)^(1-alpha)+(1-exp(u)^phi/phi)*exp(k)-exp(c);
z = mu*z(-1)+ e; 
exp(y) = exp(z)*exp(u)^alpha*exp(k)^alpha*exp(n)^(1-alpha);
end; 

% Initial values, which you should be able to figure out youselves

initval;
c=log(0.8036);
k=log(11.0836);    
n=log(0.2918);
u=log(.5);
z=1;
e=0;
y = log(1);
end;

% Variance of the shocks
shocks;
var e = 0.009^2;
end;

% Run the program
steady;
stoch_simul(periods=10100); 
