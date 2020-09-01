% Title: Performance Analysis and Non-Quadratic Lyapunov Functions for
%        Linear Time-Varying Systems
% Submitted to: American Controls Conference 2021
% Authors: Matthew Abate, Corbin Klett, Samuel Coogan and Eric Feron

% Code Author: Matthew Abate
% Date: 9/1/2020
% Description:  This script generates the Figures 2a and 2b in the paper

A = [-1 0; 0 -100];
B = [1;1];
C = [1 -2];
D = 0;
dim = 2;


dt = .001;
T = 12;
time = 0:dt:T;

y = [];
for t = time
    y = [y, expm(A*t)*b];
end

for mlevel = [1 2 5]
[Am Bm Cm] = metaSystem(A,B,C,mlevel);

cvx_begin sdp

variable P(mlevel+1,mlevel+1) semidefinite
variable gam
minimize gam
subject to
gam >= 0
[gam Cm; Cm' P] >= 0
Bm'*P*Bm <= 1
Am'*P + P*Am <= 0

cvx_end

x0m = metaState([1;1],mlevel);
x1vec = [-2:.005:2];
x2vec = x1vec;

if mlevel == 1
    P
M1 = metaPlot_2D(P, x0m,x1vec, x2vec);
bound1 = (Cm*inv(P)*Cm')^(1/(2*mlevel));
elseif mlevel == 2
M3 = metaPlot_2D(P, x0m,x1vec, x2vec);
bound3 = (Cm*inv(P)*Cm')^(1/(2*mlevel));
elseif mlevel == 5
M5 = metaPlot_2D(P, x0m,x1vec, x2vec);
bound5 = (Cm*inv(P)*Cm')^(1/(2*mlevel));
end

end


figure(1); clf; hold on; grid on; 
axis([-2 2 -1.25 1.25])
xlabel('$x_1$','Interpreter','latex');
xticks(-2:1:2)
ylabel('$x_2$','Interpreter','latex');
yticks(-1:1:1)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
set(gca, 'TickLabelInterpreter','latex','FontSize',18)


disc = 10;
% Plot impulse response ========================
plot(y(1,1:disc:end),y(2,1:disc:end), ...
            'Color', 'b', ...
            'LineWidth', 1.5);

% Plot Lyapunov function =======================
plot(M1(1,[1:disc:end, end]), M1(2, [1:disc:end, end]), ...
            'Color', 'r', ...
            'LineWidth', 1.5);

plot(M3(1,[1:disc:end, end]), M3(2,[1:disc:end, end]), ...
            'Color', [1 .5 0], ...
            'LineWidth', 1.5);
       
plot(M5(1,[1:disc:end, end]), M5(2,[1:disc:end, end]), ...
    'Color', [.2, .5, 0], ...
    'LineWidth', 1.5);
scatter(B(1), B(2), 90, 'b', 'filled')

       

% Legend ========================================
Leg = legend();
set(Leg,'visible','off')

P = [0.4099,   -0.1148; ...
   -0.1148,   0.8197]; 

arrow(1, :) = C/2.0366;
T = [cos(pi/12),sin(pi/12); -sin(pi/12), cos(pi/12)];
arrow(2, :) = C/2.0366 - .05*[1,-2]*T^1;
arrow(3, :) = C/2.0366 + .05*[1,-2]*T^11;

plot([0, C(1)/2.0366], [0, C(2)/2.0366], 'k', 'LineWidth', 1.5)
patch(arrow(:, 1), arrow(:, 2), 'k')

drawnow
%matlab2tikz('F2a.tikz', 'width', '6cm', 'height', '4cm')

%% time domain plot
h = C*y;

figure(2); clf; hold on; grid on; 
axis([0 .1 -2.6 2.6])
xlabel('$t$','Interpreter','latex');
xticks(0:.05:.1)
ylabel('$h(t)$','Interpreter','latex');
yticks(-2:1:2)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
set(gca, 'TickLabelInterpreter','latex','FontSize',18)

figure(2); hold on; grid on
plot(time(1:110),h(1:110),'Color', 'b', 'LineWidth', 1.5)
plot([0, .1],bound1*[1, 1],'Color', 'r','LineWidth', 1.5)
plot([0, .1],bound3*[1, 1],'Color', [1 .5 0], ...
            'LineWidth', 2);
plot([0, .1],bound5*[1, 1],'Color', [.2, .5, 0], ...
    'LineWidth', 2);

Leg = legend();
set(Leg,'visible','off')

plot([0, .1],-bound1*[1, 1],'Color', 'r','LineWidth', 1.5)
plot([0, .1],-bound3*[1, 1],'Color', [1 .5 0], ...
            'LineWidth', 2);
plot([0, .1],-bound5*[1, 1],'Color', [.2, .5, 0], ...
    'LineWidth', 2);

drawnow
%matlab2tikz('F2b.tikz', 'width', '6cm', 'height', '4cm')


