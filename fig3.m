% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figure 3 in the paper

clear all; clc;

% System Parameters
n = 3;
A = [.3,  0.5,   10; ...
     -1, -1.7,    1; ...
     -2,   -2, -7.7];
B = [.2; 1; 1];
C = [1 -2 2];

% Simulation Parameters
switchtime = 1; % Time at which tail is bounded
T = 2;          % Total Simulation Time
dt = .005;      % Simulation Timestep

% Simulate System and Get Impulse Response
x = B;
y = C*B;
time = 0:dt:T;
for i = 1:1:size(time, 2) - 1
    x(:, i + 1) = x(:, i)+ dt*(A*x(:, i));
    y(1, i + 1) = C*x(:, i + 1);
end

% Bound System Inpulse Response with P1
cvx_begin sdp
    variable P(n,n) semidefinite
    minimize matrix_frac(C', P);
    subject to
        B'*P*B <= 1
        A'*P + P*A <= 0
cvx_end

INITIAL_BOUND_1 = (C*inv(P)*C')^(1/2)*(B'*P*B)^(1/2);


% Bound System Inpulse Response with P4
mlevel = 4
m = n^(mlevel);
[Am Bm Cm] = metaSystem(A,B,C,mlevel);

cvx_begin sdp
    variable P(m,m) semidefinite
    minimize matrix_frac(Cm', P);
    subject to
        Bm'*P*Bm <= 10000
        Am'*P + P*Am <= 0
cvx_end

INITIAL_BOUND_4 = (Cm*inv(P)*Cm')^(1/(2*mlevel))*(Bm'*P*Bm)^(1/(2*mlevel));


% Bound The Tail of the Impulse Response Using P1
xt = x(:, switchtime/dt + 1);

cvx_begin sdp
    variable P(n,n) semidefinite
    minimize matrix_frac(C', P);
    subject to
        xt'*P*xt <= 1
        A'*P + P*A <= 0
cvx_end

TAIL_BOUND = (C*inv(P)*C')^(1/2)*(xt'*P*xt)^(1/2);


%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%

fig1 = figure(1); clf; hold on; grid on
set(gca, 'TickLabelInterpreter','latex','FontSize',18)
xlabel('$t$','Interpreter','latex');
xticks(0:1:2)
ylabel('$h(t)$','Interpreter','latex');
yticks(-1:.5:1)
axis([0 2 -1 1])

% Initial Bound From P1
plot([0, 1],  INITIAL_BOUND_1*[1, 1],'Color', 'r','LineWidth', 1.5)
plot([0, 1], -INITIAL_BOUND_1*[1, 1],'Color', 'r','LineWidth', 1.5)

% Initial Bound From P4
plot([0, 1],  INITIAL_BOUND_4*[1, 1],'Color', [.2, .5, 0],'LineWidth', 1.5)
plot([0, 1], -INITIAL_BOUND_4*[1, 1],'Color', [.2, .5, 0],'LineWidth', 1.5)

% Tail Bound From new P1
plot([1, 2],  TAIL_BOUND*[1, 1],'Color', 'r','LineWidth', 1.5)
plot([1, 2], -TAIL_BOUND*[1, 1],'Color', 'r','LineWidth', 1.5)

% Plot Imuplse Response
plot(time,y,'Color', 'b', 'LineWidth', 1.5)

Leg = legend();
set(Leg,'visible','off')

drawnow

%matlab2tikz('F3.tikz', 'width', '6cm', 'height', '4cm')
