%Macro III
%Assignment 2
%Fred Xu and Jonah Coste
%1/24/19

%Question 1
%Uses steady function defined in steady.m file.
guess = [1;1];
 
%prints answer
 steady_state = fsolve('steady', guess)
 
 %Question 2
 %2.1
 % Define E as a row matrix with two variables.
 % Each variable is an indicator for employed or unemplyed i.e.
 % E=[1,0] if employed and E=[0,1] if unemployed.
 % This is the state variable.
 % Value function becomes Value(E) = EV where V is a 1x2 column vector:
 % V = [value if employed ; value if unemployed]
 % EV = EU + beta*ETV
 % Where T is a transition matrix:
 % T = [pi_e, 1-pi_e; 1-pi_u, pi_u]
 % and U is the column matrix U=[u(w_h); u(w_l)]
 %2.2
 %Parameters
 w_h = 1;
 w_l = .9;
 beta = .96;
 pi_e = .95;
 pi_u = .9;
 
 T = [pi_e, 1-pi_e; 1-pi_u, pi_u];
 U = [log(1); log(.9)];
 I = eye(2);
 
 V=inv(I-beta*T)*U;
 
 %print answer
 V
 
 %2.3
 %Uses steady function defined in steady2.m file.
 %CAUTION beta, pi_e, and pi_u redifined in steady2.m
 guess2 = .2;
 %prints answer
 unemployment_rate = fsolve('steady2', guess2)
 
 
 %Question 3
 %nstar(z) = (w/(z*alpha))^(1/(alpha-1)
 
 % Let S_it = [s_it1; s_it2; ... ; s_itN] be a Nx1 column vector
 % where s_itx=1 if z_it = z_x and 0 otherwise.
 % Therefore S_it*Z = z_it
 % S_it is the state variable.
 % The value function is: 
 %Value(S) = SV = Spi + beta(1-lambda)STV
 % Where T is the NxN transition matrix ie
 % T_jk = f(z_i,t+1 = k | z_it = j)
 % and pi is the optimal profit row vector [pi_1, pi_2, ... , pi_N]
 % where pi_x = z_x*nstar(z_x)^alpha - w*nstar(z_x)
 
 %Question 4
 %4.1
 %Parameters
 N = 5;
 p = .8;
 w = 1;
 alpha = .7;
 beta = .95;
 lambda = .1;
 
 Z = zeros(1,N);
 T = zeros(N);
 nstar = zeros(N,1);
 pi = zeros(N,1);
 I = eye(N);
 T=T+(1-p)/(N-1)+eye(N)*(p-(1-p)/(N-1));
 for i =1:N
     Z(i) = i/N;
     nstar(i) = (w/(Z(i)*alpha))^(1/(alpha-1));
     pi(i)=Z(i)*nstar(i)^alpha - w*nstar(i);
 end
 
 %prints answer
 V = inv(I-beta*(1-lambda)*T)*pi
 
 %4.2
 Vguess = zeros(5,1);
 for i= 1:1000
     Vguess = pi + beta*(1-lambda)*T*Vguess;
 end
 
 %prints final interation
 Vguess
 
 %4.3
 iterations = 100;
 loss = zeros(iterations,1);
 Vguess = zeros(5,1);
 for j =1:5
     loss(1) = loss(1) + abs(V(j)-Vguess(j));
 end
 for i= 1:iterations-1
     Vguess = pi + beta*(1-lambda)*T*Vguess;
     for j =1:5
        loss(i+1) = loss(i+1) + abs(V(j)-Vguess(j));
     end
 end
 
 plot(loss)
 title('Value Function Convergence--True Value Known')
 xlabel('Iteration')
 ylabel('Loss (sum of absolute differences from true value')
 
 % ~67 iterations gives a fairly precise estimate of the value function.
 % Note this is the first iteration where loss function rounds to .0000
 %4.4
 loss2 = zeros(iterations,2);
 Vold = zeros(5,1);
 Vnew = zeros(5,1);
 for i= 1:iterations
     loss2(i,1) = loss(i);
     Vold=Vnew;
     Vnew = pi + beta*(1-lambda)*T*Vold;
     for j =1:5
        loss2(i,2) = loss2(i,2) + abs(Vnew(j)-Vold(j));
     end
 end
 
 plot(loss2)
 title('Value Function Convergence')
 xlabel('Iteration')
 ylabel('Loss (sum of absolute differences')
 legend({'True Value Known','True Value Unknown'},'Location','southeast')
 
 % Now ~55 iterations gives a fairly precise estimate of the value function.
 % Note this is the first iteration where loss function rounds to .0000
