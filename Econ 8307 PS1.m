%Macro III
%Assignment 1
%Fred Xu and Jonah Coste
%1/17/19

%Question 1
%1.A
%1.A.1
N = randi(10);
%1.A.2
X = rand(N,1)*2 -1;
%Expected value of mean is 0.
%1.A.3.1
sum1 = 0;
for i=1:N
    sum1 = sum1 + X(i);
end
result1 = sum1/N;
%1.A.3.2
 result2 = mean(X);
 %1.A.4
 result1
 result2
 
 %1.B
 %1.B.1
 N = randi(10);
 M = randi(10);
 L = randi(10);
 %1.B.2
 matrix1 = rand(N,M)*2 -1;
 matrix2 = rand(M,L)*2 -1;
 %1.B.3.1
 resultmatrix1 = zeros(N,L);
 for i =1:N
     for j = 1:L  
        val = 0;
        for h = 1:M
            val = val + matrix1(i,h) * matrix2(h,j);
        end
        resultmatrix1(i,j) = val;
     end
 end
 %1.B.3.2
 resultmatrix2 = matrix1 * matrix2;
 %1.B.4
 resultmatrix1
 resultmatrix2
 
 %Question 2
 %2.A
 %Use World Bank Indicator data for GDP in constant 2010 US$ for China and
 %USA for years 1960-2017 available here https://data.worldbank.org/indicator/NY.GDP.PCAP.KD?locations=CN-US
 %2.B
 China = [191.791179910216;140.913499235072;131.96337147674;142.021955444481;163.991325226585;187.274631568573;201.52324200402;185.075927473726;172.913752439848;196.740126889365;228.317702823625;237.813838393538;240.881888917967;253.714372960062;254.267484665894;271.599476492383;263.23062210942;279.32454675135;307.76619470515;326.768369193998;347.887413035502;361.224710649125;387.745580908334;423.593498767967;481.364595759347;538.690827266253;578.184040366224;635.494603029004;695.599054411093;713.689527590362;730.772489044822;787.867435156269;888.911004119008;1000.61180975351;1118.49957748173;1227.55640691521;1335.36268011202;1443.77474185404;1542.06412996664;1645.98799567879;1771.74150579539;1905.61078010947;2065.71857925612;2258.91210541049;2472.58655569402;2738.2054599526;3069.30478095137;3487.84576561013;3805.02599866378;4142.03828597868;4560.51258600929;4971.544928635;5336.06014319861;5721.69381888796;6108.23877494863;6496.62401255517;6894.46452231334;7329.08929913216];
 USA = [17036.8851695882;17142.1937673897;17910.2787901684;18431.1584040472;19231.1718590604;20207.7495376936;21274.1354890107;21569.8356088859;22380.6067673191;22850.010833783;23309.6209459064;23775.2769229658;24760.1453772345;25908.9128017215;25540.5010030208;25239.919905729;26347.809281996;27286.2515144911;28500.2404573534;29082.5937779654;28734.3992597646;29191.9994879416;28362.4946163409;29406.2574686046;31268.9756446518;32306.8330567744;33133.6954441913;33975.6547953064;35083.9690428182;36033.3302026992;36312.4141825874;35803.8684213439;36566.1737698527;37078.04968394;38104.9724675631;38677.7150883663;39681.5198579033;40965.8466450522;42292.8912011426;43768.8849928326;45055.817918284;45047.4871976844;45428.6456781274;46304.0360895612;47614.2798621765;48755.6160606735;49575.401013591;49979.5338429195;49364.6445500336;47575.608562749;48375.4069462972;48786.4549755253;49498.3909155226;49971.9513571108;50871.674083306;51933.4048064982;52319.1633505427;53128.5396999252];
 %2.C
 sz = size(China,1) - 1;
 Growth = zeros(sz ,2);
 Year = zeros(sz,1);
 for i = 1:sz
     Growth(i,1) = China(i+1,1) / China(i,1) -1;
     Growth(i,2) = USA(i+1,1) / USA(i,1) -1;
     Year(i,1) = 1960 + i;
 end
 figure
 plot(Year, Growth)
 title('Real GDP Growth Rate for China and USA between 1960-2017')
xlabel('Year')
ylabel('Growth Rate')
legend({'China','USA'},'Location','southeast')
 
 %Question 3
 %3.A
 % Agent's budget constraint is:
 % f(k_t) + (1-delta) * k_t >= c_t + k_{t+1}
 % aslo k_{t+1} >= 0
 %3.B
 % From t=1 to T-1 the first order condition that charachterizes the
 % solution is:
 % u'(c_t) = beta * u'(c_{t+1}) * [f'(k_{t+1}) + 1 - delta]
 % The optimality condtion at period T is:
 % c_T =  f(k_T) + (1-delta) * k_T
 % At the end of period T there should be zero capital left over (i.e.
 % k_{T+1} = 0)
 % 3.C
 %parameters
 beta = .97;
 delta = .1;
 theta = .3;
 T = 4;
 epsilon=1e-15;
 initialcapital = 1;
 guessinitialcons = 0; %could be any guess
 
 if guessinitialcons < 0
     guessinitialcons = 0;
 end %fail safe
 
 c = zeros(T,1);
 k = zeros (T+1,1);
 k(1) = initialcapital;
 
 guess = guessinitialcons;
 last = -2*epsilon; % to allow for guess of 0
 step = k(1);
 
 while abs(guess - last) > epsilon || abs(k(T+1)) > epsilon
     c(1) = guess;
     ind = 0;
     for i= 1: T-1
        if k(i)^theta + (1 - delta)*k(i) - c(i) >= 0
            k(i+1) = k(i)^theta + (1 - delta)*k(i) - c(i);
            c(i+1) = beta * c(i) * (theta * k(i+1)^(theta-1) + 1 - delta);
        else
            ind = 1;
        end
     end
     k(T+1) = k(T)^theta + (1 - delta)*k(T) - c(T);
     if k(T+1) < 0 || ind == 1
         step = step/2;
         guess = last + step;
     else
         last = guess;
         guess = guess + step;
     end
 end
 
 %Print optimal consumption and capital series
 c
 k
 
 
