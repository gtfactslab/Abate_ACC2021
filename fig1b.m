% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figure 1b in the paper
clc; clear all;

A = [0, 1; -1, -0.9];
b = [1; 1];
c = [sqrt(2), -sqrt(2)];

dt = .01;
T = 12;
time = 0:dt:T;

impulse = [];
for t = time
    impulse = [impulse, c*expm(A*t)*b];
end


% figure 1
figure(1); clf; hold on; grid on; 
axis([0 T -3 3])
xlabel('$t$','Interpreter','latex')
xticks(0:T/4:T)
ylabel('$h(t)$','Interpreter','latex')
yticks(-3:1:3)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

plot([0, T], 2*sqrt(2)*[1, 1], 'r', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([0, T], -2*sqrt(2)*[1, 1], 'r', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(time(1:10:end), impulse(1, 1:10:end), 'b', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');

Leg = legend();
set(Leg,'visible','off')
drawnow

% matlab2tikz('F1b.tikz', 'width', '6cm', 'height', '4cm')
