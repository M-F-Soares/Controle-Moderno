clear all; close all; format short; clc;

m1 = 1; m2 = 1; l1 = 1; l2 = 1; g = 9.81;

x = sym('x', [1 4], 'real');

f1 = ((sin(x(1))*cos(x(3)) - cos(x(1))*sin(x(3)))*(m2*l1*l2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))*(x(2)^2)-m2*(l2^2)*(x(4)^2)))/(l1*l2*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2)) + ((m1+m2)*l2*g*sin(x(1))-m2*l2*g*sin(x(3))*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3))))/(l1*l2*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2));
f2 = ((sin(x(1))*cos(x(3)) - cos(x(1))*sin(x(3)))*(-(m1+m2)*(l1^2)*(x(2)^2)+m2*l1*l2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))*(x(4)^2)))/(l1*l2*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2)) + (-(m1+m2)*l1*g*sin(x(1))*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))+(m1+m2)*l1*g*sin(x(3)))/(l1*l2*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2));

g11 = (m2*l2^2)/(m2*(l1^2)*(l2^2)*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2));
g12 = (-m2*l1*l2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3))))/(m2*(l1^2)*(l2^2)*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2));
g21 = g12;
g22 = ((m1+m2)*(l1^2))/(m2*(l1^2)*(l2^2)*((m1+m2)-m2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)))^2));

dotx = sym('dx', [1 4]);
tau = sym('tau', [1 2], 'real');

dotx(1) = x(2);
dotx(2) = f1 + g11*tau(1) + g12*tau(2);
dotx(3) = x(4);
dotx(4) = f2 + g21*tau(1) + g22*tau(2);

y = sym('y', [1 2], 'real');
y(1) = x(1); %theta 1
y(2) = x(3); %theta 2

%Definição dos ângulos que definem o ponto de operação (yop)
X_op(1) = pi/2; %ângulo da junta 1
X_op(2) = 0;
X_op(3) = -pi/2; %ângulo da junta 2
X_op(4) = 0;
for i=1:4
    X(i) = subs(dotx(i), {x(1), x(2), x(3), x(4)}, {X_op(1), X_op(2), X_op(3), X_op(4)});
end

%Cálculo do sinal de controle do ponto de operação (uop)
aux1 = vpa(solve(subs(X(2)==0,{x(2) x(4)}, {X_op(2) X_op(4)}), tau(1)), 5);
Tau_op(2) = vpa(solve(subs(X(4), {x(2) x(4) tau(1)}, {X_op(2) X_op(4) aux1})==0, tau(2)), 5); %torque da junta 2
Tau_op(1) = vpa(subs(aux1, {x(2) x(4) tau(2)}, {X_op(2) X_op(4) Tau_op(2)}), 5); %torque da junta 1

for i=1:size(dotx,2)
    for j=1:size(dotx,2)
        if j == 1
            B(i, 1) = subs(diff(dotx(i), tau(1)), {tau(1), tau(2), x(1), x(2), x(3), x(4)},{Tau_op(1), Tau_op(2), X_op(1), X_op(2), X_op(3), X_op(4)});
        else
            B(i, 2) = subs(diff(dotx(i), tau(2)), {tau(1), tau(2), x(1), x(2), x(3), x(4)},{Tau_op(1), Tau_op(2), X_op(1), X_op(2), X_op(3), X_op(4)});
        end
        A(i, j) = subs(diff(dotx(i), x(j)),{tau(1), tau(2), x(1), x(2), x(3), x(4)},{Tau_op(1), Tau_op(2), X_op(1), X_op(2), X_op(3), X_op(4)});
    end     
end

A = double(vpa(A, 5));
B = double(vpa(B, 5));
C = [1 0 0 0; 0 0 1 0];
D = zeros(size(C,1), size(C,1));
x = x';
u = tau';
Tau_op = double(Tau_op)';

sis = ss(A,B,C,D)
step(sis)
grid

%%
% ESPECIFICAÇÔES DE DESEMPENHO %
OS = 5;
ts = 1;
zeta = -log(OS/100)/(sqrt(pi^2 + (log(OS/100))^2));
wn = 4/(zeta*ts);
wbw = wn*sqrt((1-2*zeta^2) + sqrt(4*zeta^4+4*zeta^2+2));

Dom_poles = [-zeta*wn + wn*sqrt(zeta^2-1) -zeta*wn - wn*sqrt(zeta^2-1)];
for i=1:(size(A,1)+ size(D,1))
    if i==1 
        P(i) = Dom_poles(1);
    elseif i==2
        P(i) = Dom_poles(2);
    else
        P(i) = real(Dom_poles(1)*10-(i-2));
    end
end
%n+1 pólos graças à presinça da ação integral na matriz Aa

%%
% REALIMENTAÇÃO DE ESTADOS COM AÇÃO INTEGRAL %
%Confere se o sistema é controlável
if(size(A,1)==rank(ctrb(A,B)))
    disp('Controlável')
end

Aa = [A zeros(size(A,1),2);-C zeros(size(D))]; %A aumentada pela ação integral
Ba = [B; zeros(size(D))]; %B aumentada pela ação integral

K_ = [1 0 1 0 1 0; 0 1 0 1 0 1];

for i=1:(size(A,1)+ size(D,1))
    for j=1:(size(A,2)+ size(D,2))
        if (i==1 && j==1) || (i==2 && j==2)
            Fctr(i,j) = real(Dom_poles(i));
        elseif (i==2 && j==1) || (i==1 && j==2)
            Fctr(i,j) = imag(Dom_poles(i));
        elseif i==j
            Fctr(i,j) = P(i);
        else
            Fctr(i,j) = 0;
        end
    end
end

%Confere se o par (Fctr, K_) é observável
if(size(Fctr,1)==rank(obsv(Fctr,K_)))
    disp('Observável')
end

Tctr = lyap(Aa, -Fctr, -Ba*K_);
Kt = K_/Tctr;
K = [Kt(1, 1:size(A,1)); Kt(2, 1:size(A,1))];
Ka = [-Kt(1,end-(size(B,2)-1):end); -Kt(2, end-(size(B,2)-1):end)]; 
    
%%
% OBSERVADOR %
if(size(A,1)==rank(obsv(A,C)))
    disp('Observável')
end

L = [1 0; 0 1; 1 0; 0 1];

for i=1:(size(A,1))
    for j=1:(size(A,2))
        if i==j
            Fobs(i,j) = real(P(1))*5-(i-1);
        else
            Fobs(i,j) = 0;
        end
    end
end

%Confere se o par (Fobs,L) é controlável
if(size(Fobs,1)==rank(ctrb(Fobs,L)))
    disp('Controlável')
end

Tobs = lyap(-Fobs, A, -L*C);