% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figure 8a in the paper

clear
close all;

% This file to compute the max difference between two outputs
Anom = [0 1;-0.6 -0.5];
Ad = [0 0;0.1 -0.1];
A11 = Anom+Ad;
A22 = Anom-Ad;

B = [0;1];
C = [1 0];
[n,o]=size(A11);
n = 2*n;

% Compute bounds
A1 = [A11 zeros(2,2) ; zeros(2,2) Anom];
A2 = [A22 zeros(2,2) ; zeros(2,2) Anom];
C1 = [C -C];
B1 = [B;B];

cmax = 3; %orders of MLFs (DO NOT INCLUDE A VALUE OVER 6)

% Largest
div = [0.1];

%%% Plot MLF method
BOUND = zeros(cmax,1);
for c = 1:cmax %order of Lyapu = 2c
   [P,gammma] = commonLyap(c,A1,A2,B1,C1,n, div);
   BOUND(c)=gammma^(1/(2*c));
end

BOUND

%
T = 12;
dt = .01;
time = 0:dt:T;

figure(1); clf; hold on; grid on
set(gca, 'TickLabelInterpreter','latex','FontSize',18)
xlabel('$t$','Interpreter','latex');
xticks(0:T/2:T)
ylabel('$h(t)$','Interpreter','latex');
yticks(-.75:.75:1.5)
axis([0 T -.75 1.5])

%some cheap simulations
Q = {Anom, A11, A22};
for i = 1:3
    [num,den] = ss2tf(Q{1, i},B,C,0);
    y=impulse(num,den,time);
    plot(time,y, 'b', 'LineWidth', 1.5);
end

[num,den] = ss2tf(Anom,B,C,0);
y=impulse(num,den,time);

%%%%%%%%%%%%
% Plot Figure 8a:
%%%%%%%%%%%%

% Plot bound alpha_1_max at P1 level
boundu = y + BOUND(1)*exp(-div*time)';
boundl = y - BOUND(1)*exp(-div*time)';
plot(time, boundu, 'Color', [1,  0, 0], 'LineWidth', 1.5);
plot(time, boundl, 'Color', [1,  0, 0], 'LineWidth', 1.5);

% Plot bound alpha_2_max at P2 level
boundu = y + BOUND(2)*exp(-div*time)';
boundl = y - BOUND(2)*exp(-div*time)';
plot(time, boundu, 'Color', [1, .5, 0], 'LineWidth', 1.5);
plot(time, boundl, 'Color', [1, .5, 0], 'LineWidth', 1.5);

% Plot bound alpha_2_max at P2 level
boundu = y + BOUND(3)*exp(-div*time)';
boundl = y - BOUND(3)*exp(-div*time)';
plot(time, boundu, 'Color', [.2,  .5, 0], 'LineWidth', 1.5);
plot(time, boundl, 'Color', [.2,  .5, 0], 'LineWidth', 1.5);

Leg = legend();
set(Leg,'visible','off')

drawnow

%matlab2tikz('F8a.tikz', 'width', '6cm', 'height', '4cm')


function [out,out1] = commonLyap(c, A1,A2, B, C,n,div)
A1 = A1+div*eye(4);
A2 = A2 +div*eye(4);
if c == 1
        Ac1 = A1
        Ac2 = A2
        Bc = B;
        Cc=C;

    else
        Ac1 = 0;  %A2c = 0;
        Ac2 = 0
        Bc = B;
        Cc = C;
        for i = 0:c - 1
            Ac1 = Ac1 + kron(kron(eye(n^(i)), A1), eye(n^(c-1 - i)));
            Ac2 = Ac2 + kron(kron(eye(n^(i)), A2), eye(n^(c-1 - i)));
            if i>0 
                Bc = kron(Bc,B);
                Cc = kron(Cc,C);
            end
        end
    end

    %Search for common Lyapunov function
 
    cvx_begin sdp

        %Optimize over P
        variable P(n^(c),n^(c)) semidefinite
        variable gammma

        %Constraints
        0 >= Ac1'*P + P*Ac1;
        0 >= Ac2'*P + P*Ac2;
        
        Bc'*P*Bc <= 1;
        [P Cc';Cc gammma] >= 0;
        
        minimize(gammma);
        
    cvx_end;
    if isequal(cvx_status, 'Failed') || isequal(cvx_status, 'Infeasible')
        gammma=0;
    end
    out = P;
    out1 = gammma;
 
end
