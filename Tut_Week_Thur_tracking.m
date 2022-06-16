%% Ex,3.2-3 (p.118) - free final state
% J(0)  = 1/2 s(x(T) - 10)^2 ) + 1/2 Integral _0^T (u^2 )dt
clear all; clc; clf
LW='Linewidth';
T = 10;
a =1; b=1;
S= [1 10, 100];

for i  =1:3
    x = chebfun(@(t) 10*S(i)*b^2*sinh(a*t) /(a*exp(a*T) + S(i)*b^2*sinh(a*T)), [0, T]);
    u = chebfun(@(t) 10*a*b*S(i)*exp(a*t)/(a*exp(a*T) + S(i)*b^2*sinh (a*T)), [0, T]);
    figure(1)
    plot(u,LW,2); grid on; hold on
    title(' optimal controller  with S =[1 10 100]')
    figure(2)
    plot(x,LW,2); grid on ; hold on
    title(' optimal temperature  with S =[1 10 100]')
end 
hold off


%% Ex.3.2-3, chebfun_guide: ch.7.8

clear all; clc;
LW='Linewidth';
T = 10;
a =1; b=1;
S=100;
L = chebop(0,10);
L.op=@(x,u,v) [diff(u)+u+v; diff(v)-v];
L.lbc =@(u,v) u;
L.rbc =@(u,v) u*S-v-10*S;
rhs =[0;0];
U =L\rhs;
figure(3)
plot(U(1,:),LW,2); grid on


%% Ex3.3-5
clear all;clc;clf
om=0.8;    del=0.1;   
a =[0 1; -om^2 -2*del*om];
b=1;r=1; x0 = [10 0]';
[x,u,Sf,tf]=ex3_3_5(a,b,r,x0);
% whos

figure(1)
plot(tf,x(1,:)); grid on;
figure(2)
plot(tf,u); grid on

%% Ex. 4.1-1 - p.180
clear all; clc;clf
LW = 'linewidth',
a = -1; b=1; T = 100;
p = 10;  R = 1;
D =[0, T];
h =1;
r = chebfun(@(t) h*sign(t), D,'splitting', 'on');
% r = chebfun(@(t) sin(t), D,'splitting', 'on');
q1=[10 100 1000];

for i = 1:3
  q=q1(i);
    S = chebop(D); S.rbc =p;
    S.op = @(t,s) diff(s)+2*a*s-b^2*s^2/R + q;
    s =S\0; 
    figure(1) ;
    plot(s, LW,2); grid on ; hold on
    V = chebop(D); V.rbc =p*r(T);
    V.op=@(t,v) diff(v) +(a-b^2*s/R)*v + q*r;
    v = V\0;
    figure(2)
    plot(v,LW,2); grid on; hold on
    K = b*s/R;
    X = chebop(D); X.lbc =0;
    X.op = @(t,x) diff(x) -a*x + b*K*x - b*v/R;
    x = X\0;
    figure(3)
    plot(x,LW,2); grid on; hold on
    plot(r,LW,2,'r')
   
    % cost evaluation
    % the final value of w(t)
    W_T = 1/2*p*(x(T)-h)^2-1/2*x(T)^2*p + x(T)*v(T);
    W = chebop(D); W.rbc =W_T;
    W.op=@(t,w) diff(w)+1/2*r^2*q-1/2*v^2*b^2/R;
    w=W\0;
    J = 1/2*x^2*s-x*v+w;
    w=W\0;
    figure(4)
    plot(J,LW,2); grid on; hold on 
    pause(2)
end 
hold off

%% Numerical solution 
clear all; clc;
L=chebop(0, 10);
L.op=@(x,u,v) [diff(v)+u; diff(u)-v];
L.lbc =@(u,v) u;
L.rbc =@(u,v) v;
rhs =[1,1]';
U=L\rhs;
plot(U)



%% Ex. 4.1-2 p.182 : disturbance rejection
% dx = ax + bu + d
% J = 1/2 x^2 Sx^2 + 1/2 integral(x^TQx + u^TRu)
% boundary condition
% BrySon , P.176

clear all; clc;clf
LW = 'linewidth',
a = -1; b=1; T = 10;
p = 10;  R = 1;
D =[0, T];
%r = chebfun(@(t) sign(t), D,'splitting', 'on');
d = chebfun(@(t) sin(t), D);
q1=[10 100 1000];

for i = 1:3
  q=q1(i);
    S = chebop(D); S.rbc =p;
    S.op = @(t,s) diff(s)+2*a*s-b^2*s^2/R + q;
    s =S\0; 
    figure(1) ;
    plot(s, LW,2); grid on ; hold on
    K = b*s/R;
    V = chebop(D);  V.rbc =0;
    V.op=@(t,v) diff(v) +(a-b*K)*v + s*d;
    v = V\0;
    figure(2)
    plot(v,LW,2); grid on; hold on
    X = chebop(D); X.lbc =0;
    X.op = @(t,x) diff(x) -a*x + b*K*x + b^2*v/R;
    x = X\0;
    figure(3)

    plot(x,LW,2); grid on; hold on
    plot(d,LW,2,'r')
    pause(2)
    % plot(-K*x+b*v/R, LW,2,'+b')
  
end
hold off
%% Kim. Ex4.1-3 final state fixed
clear all; clc;clf
LW = 'linewidth';
a = -1; b=1; T = 10; c= 1; r(T) = 0;
p = 10;  R = 1;
D =[0, T];
%r = chebfun(@(t) sign(t), D,'splitting', 'on');
% d = chebfun(@(t) sin(t), D,'splitting', 'on');
q1=[10 100 1000];

for i = 1:1
  q=q1(i);
    S = chebop(D); S.rbc =p;
    S.op = @(t,s) diff(s)+2*a*s-b^2*s^2/R + q;
    s =S\0; 
    figure(1) ;
    plot(s, LW,2); grid on ; hold on
    K = b*s/R;
    V = chebop(D);  V.rbc =c;
    V.op=@(t,v) diff(v) +(a-b*K)*v ;
    v = V\0;
    figure(2)
    plot(v,LW,2); grid on; hold on
    P = chebop(D); P.lbc =0;
    P.op =@(t,p) diff(p)-(b^2*v^2)/R;
    p=P\0;
    figure(3)
    plot(p,LW,2); grid on; hold on
%     X = chebop(D); X.lbc =0;
%     X.op = @(t,x) diff(x) -a*x + b*v^2/(R*p)*x + b*v*r(T)/(R*p);
%     X.op = @(t,x) diff(x) -a*x  + b*v*r(T)/(R*p);
%     x = X\0;
%     figure(4)
%     plot(x,LW,2); grid on; hold on
%     plot(d,LW,2,'r')
%     pause(2)
    % plot(-K*x+b*v/R, LW,2,'+b')
  
end
hold off



%%  Exploring.ODE: systems of equations : page 116
clear all; clc;
LW = 'linewidth',
a = -1; b=1; T = 10; c= 1; r(T) = 0;
p = 10;  R = 1;q=1;
D =[0, T];
r =0;
N = chebop(D); N.lbc=[1;0];
N.op = @(t,x,y) [diff(x) - a*x +(b/R)*y; diff(y) +a* y + c^2*q*x];
N.lbc=[1;0];
[x,y] = N\0;
plot([x,y]); grid on;
    



