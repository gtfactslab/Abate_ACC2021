
% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figure 1a in the paper


clc; clear all;

A = [0, 1; -1, -0.9];
b = [1; 1];
c = [sqrt(2), -sqrt(2)];

impulse = [];
dt = .01;
time = 0:dt:12;

holder = 0;
for t = time
    x_now = expm(A*t)*b;
    impulse = [impulse, x_now];
    
    if abs(c*x_now) >= abs(c*holder)
        holder = x_now
    end
end
impulse_responce = holder;



 %Search for common Lyapunov function using Algorithm 1
cvx_begin sdp
    %Optimize over P
    variable P(2, 2) semidefinite
    
    %Constraints
    1 >= b.'*P*b;
    0 >= A.'*P + P*A;
    minimize(matrix_frac(c', P));
cvx_end
P

% figure 1
figure(1); clf; hold on; grid on; 
axis([-2 2 -1.5 1.5])
xlabel('$x_1$','Interpreter','latex')
xticks(-1.5:.5:1.5)
ylabel('$x_2$','Interpreter','latex')
yticks(-1:.5:1)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


plot(sqrt(2)*sin(0:.01:2*pi), sqrt(2)*cos(0:.01:2*pi), 'r', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');

plot(holder(1) + (-1:.1:1), holder(2)+ (-1:.1:1), 'm--')
plot(1 + (-1:.1:1), -1+ (-1:.1:1), 'm--')

scatter(b(1), b(2), 90, 'b', 'filled')
plot(impulse(1, :), impulse(2, :), 'b', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');

plot([0, c(1)], [0, c(2)], 'k', 'LineWidth', 1.5)

scatter(1, -1, 90, 'r', 'filled')

scatter(impulse_responce(1), impulse_responce(2), 90, 'm', 'filled')

patch(c(1) + [-.15, 0, -.1] +.03, c(2) + [.1, 0, .15] -.03, 'k')

Leg = legend();
set(Leg,'visible','off')
drawnow

% matlab2tikz('F1a.tikz', 'width', '7cm', 'height', '5cm')
