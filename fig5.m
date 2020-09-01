% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figure 5 in the paper

close all; clear;

Anom = [0 1;-0.6 -0.5];
Ad = [0 0;0.1 -0.1];
div=[0.0;0.16;-0.5];
A1 = Anom + Ad;
A2 = Anom - Ad;

B = [0;1];
C = [1 0];

[n,o]=size(Anom);

T = 10;
dt = .01;
time = 0:dt:T;



% Compute bounds
cmax = 6; %orders of MLFs (DO NOT INCLUDE A VALUE OVER 6)

%%% Plot MLF method
BOUND = zeros(cmax,2);
for c = 1:cmax %order of Lyapu = 2c
    [Am1 Am2 Bm Cm] = metaSystem2(A1, A2, B,C,c);
      
    cvx_begin sdp

        %Optimize over P
        variable P(n^(c),n^(c)) semidefinite

        %Constraints
        0 >= Am1'*P + P*Am1;
        0 >= Am2'*P + P*Am2;

        Bm'*P*Bm <= 1;
        minimize(matrix_frac(Cm', P));
        
    cvx_end
    Q{1, c} = P;
    BOUND(c, 1)=(Cm*inv(P)*Cm')^(1/(2*c))*(Bm'*P*Bm)^(1/(2*c)); 
end



alpha = -.5;
A1a = A1 + alpha*eye(2);
A2a = A2 + alpha*eye(2)

for c = 1:cmax %order of Lyapu = 2c
    [Am1 Am2 Bm Cm] = metaSystem2(A1a, A2a, B,C,c);
      
    cvx_begin sdp

        %Optimize over P
        variable P(n^(c),n^(c)) semidefinite

        %Constraints
        0 >= Am1'*P + P*Am1;
        0 >= Am2'*P + P*Am2;

        Bm'*P*Bm <= 1;
        minimize(matrix_frac(Cm', P));
        
    cvx_end
    Q{2, c} = P;
    BOUND(c, 2)=(Cm*inv(P)*Cm')^(1/(2*c))*(Bm'*P*Bm)^(1/(2*c)); 
end


beta = .15;
A1a = A1 + beta*eye(2);
A2a = A2 + beta*eye(2)

for c = 1:cmax %order of Lyapu = 2c
    [Am1 Am2 Bm Cm] = metaSystem2(A1a, A2a, B,C,c);
      
    cvx_begin sdp

        %Optimize over P
        variable P(n^(c),n^(c)) semidefinite

        %Constraints
        0 >= Am1'*P + P*Am1;
        0 >= Am2'*P + P*Am2;

        Bm'*P*Bm <= 1;
        minimize(matrix_frac(Cm', P));
        
    cvx_end
    Q{3, c} = P;
    BOUND(c, 3)=(Cm*inv(P)*Cm')^(1/(2*c))*(Bm'*P*Bm)^(1/(2*c)); 
end


figure(1); clf; hold on; grid on
set(gca, 'TickLabelInterpreter','latex','FontSize',18)
xlabel('$t$','Interpreter','latex');
xticks(0:T/2:T)
ylabel('$h(t)$','Interpreter','latex');
yticks(-1:1:1)
axis([0 T -1.1 1.1])




%some cheap simulations
Q = {Anom, A1, A2};
boundu = [];
boundl = [];

[num,den] = ss2tf(Anom,B,C,0);
y=step(num,den,time);

for i = 1:3
    [num,den] = ss2tf(Q{1, i},B,C,0);
    y=impulse(num,den,time);
    plot(time,y, 'b', 'LineWidth', 1.5);
end


Color = [ 1,  0, 0 ;
          1,  .5, 0; ...
          0,  0, 0; ...
          0,  0, 0; ...
          0,  0, 0; ...
          .2,  .5, 0];
      
plot([0, T],  BOUND(6, 1)*[1, 1],'Color', Color(1, :), 'LineWidth', 1.5)
plot([0, T], -BOUND(6, 1)*[1, 1],'Color', Color(1, :), 'LineWidth', 1.5)
plot(time, BOUND(end, 2)*exp(-alpha*time),'Color', Color(2, :), 'LineWidth', 1.5)
plot(time, -BOUND(end, 2)*exp(-alpha*time),'Color', Color(2, :), 'LineWidth', 1.5)
plot(time, BOUND(end, 3)*exp(-beta*time),'Color', Color(6, :), 'LineWidth', 1.5)
plot(time, -BOUND(end, 3)*exp(-beta*time),'Color', Color(6, :), 'LineWidth', 1.5)

Leg = legend();
set(Leg,'visible','off')

drawnow

%matlab2tikz('F5.tikz', 'width', '6cm', 'height', '4cm')



