% This is the file from the notes
% Remember to save files like this where MATLAB can find them

% Tell dynare what variables and parameters there are
var c, k, n, z, u, y, w; 
varexo e; 
parameters beta, alpha, theta, gam, mu, phi; 
 
% Give values to the parameters
alpha = 0.36; 
mu = 0.95; 
beta = 0.99;  
gam =0; 
theta = 2.95; 
phi = 1.2;

% The dynamic equations go here
% Note the conventions for time. State variables like capital go from t to t+1
% Endogenous variables like c and z go from t-1 to t. Labor n does not appear dynamically

model;
beta*c^(-1)*(1-u^phi/phi+alpha*exp(z)*u^alpha*k^(alpha-1)*n^(1-alpha))=c(-1)^(-1);
theta*(1-n)^(-gam)=c^(-1)*(1-alpha)*exp(z)*u^alpha*k^alpha*n^(-alpha);
alpha*exp(z)*u^(alpha-phi)*k^(alpha-1)*n^(1-alpha)=1;
k(+1)=exp(z)*u^alpha*k^alpha*n^(1-alpha)+(1-u^phi/phi)*k-c;
z = mu*z(-1)+ e; 
y = exp(z)*u^alpha*k^alpha*n^(1-alpha);
w = mu*z(-1)+ f;
end; 

% Initial values, which you should be able to figure out youselves

initval;
c=0.8036;
k=11.0836;    
n=0.2918;
u=.5;
z=0;
e=0;
y = 1;
end;

% Variance of the shocks
shocks;
var e = 0.009^2;
end;

% Run the program
steady;
stoch_simul(periods=10100); 
