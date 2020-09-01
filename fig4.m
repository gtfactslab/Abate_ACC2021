% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figure 4 in the paper

clc; clear all;

A = [-1 0 2; 0 -10 1; 0 -2 -1];
B = [-2;1;1];
C = [1 -2 2];


n = size(A,1);

T = 6;
switch_time = 2;
dt = .005;
x = [0; 0; 0];
y = 0;

time = 0:dt:T;
for i = 1:1:size(time, 2) - 1
    x(:, i + 1) = x(:, i)+ dt*(A*x(:, i) + B);
    y(1, i + 1) = C*x(:, i);
end


LEVELS = 1:3;
for mlevel = LEVELS
    m = n^(mlevel);
    [Am Bm Cm] = metaSystem(A, -inv(A)*B,C,mlevel);
    bound = step_bound(Am, Bm, Cm)

    b(mlevel, 1) = (bound)^(1/(2*mlevel));
end


for mlevel = LEVELS
    m = n^(mlevel);

    xt = x(:, switch_time/dt+1);
    [Am Bm Cm] = metaSystem(A, -xt-inv(A)*B,C,mlevel);

    bound = step_bound(Am, Bm, Cm)
    b(mlevel, 2) = (bound)^(1/(2*mlevel));
end

q = -C*inv(A)*B;
%%
fig1 = figure(1); clf; hold on; grid on
set(gca, 'TickLabelInterpreter','latex','FontSize',18)
xlabel('$t$','Interpreter','latex');
xticks(0:2:6)
ylabel('$h(t)$','Interpreter','latex');
yticks(-.5:.5:1)
axis([0 T -.5 1.25])

Color = [ 1,  0, 0 ;
          1,  .5, 0; ...
          .2,  .5, 0];
for i = [1, 3]
    plot([0, switch_time],  b(i, 1)*[1, 1] + q,'Color', Color(i, :), 'LineWidth', 2)
    plot([0, switch_time], -b(i, 1)*[1, 1] + q,'Color', Color(i, :), 'LineWidth', 2)
end
for i = 1:1
    plot([switch_time, T],  b(i, 2)*[1, 1] + q,'Color', Color(i, :), 'LineWidth', 2)
    plot([switch_time, T], -b(i, 2)*[1, 1] + q,'Color', Color(i, :), 'LineWidth', 2)
end

plot(time,y,'Color', 'b', 'LineWidth', 1.5)

Leg = legend();
set(Leg,'visible','off')

drawnow

%matlab2tikz('F4.tikz', 'width', '7cm', 'height', '5cm')
