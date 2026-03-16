%%ECE 141 Final Project -- Wesley Hon

%Make sure you home everything.
clear; clc; close all;

%% PART 0: Parameters
R   = 1000;          % ohms
P   = 20;            % watts
L   = 100e-6;        % henries
C   = 100e-6;        % farads
vs0 = 80;            % equilibrium source voltage (V)

%% PART 1: Nonlinear state-space representation
syms x1 x2 u real

f1 = (1/L)*(u - x2);
f2 = (1/C)*(x1 - x2/R - P/x2);
f  = [f1; f2];

fprintf('Part 1: Nonlinear State-Space Model\n');
disp('f(x,u) = ');
pretty(f)

%% PART 2: Linearization
% Jacobians
x = [x1; x2];
A_sym = jacobian(f, x);
B_sym = jacobian(f, u);

fprintf('\nPart 2: Symbolic Linearization\n');
disp('A(x,u) = ');
pretty(A_sym)
disp('B(x,u) = ');
pretty(B_sym)

%% PART 3 & 4: Equilibrium point
vdc0 = vs0;
iL0  = vdc0/R + P/vdc0;

fprintf('\nPart 4: Equilibrium Values\n');
fprintf('vdc0 = %.4f V\n', vdc0);
fprintf('iL0  = %.4f A\n', iL0);

%% PART 2 continued: Evaluate A and B at equilibrium
A = double(subs(A_sym, [x1 x2 u], [iL0 vdc0 vs0]));
B = double(subs(B_sym, [x1 x2 u], [iL0 vdc0 vs0]));

Cout = [0 1];   % output y = delta x2 = delta vdc
D = 0;

fprintf('\nLinearized Matrices at Equilibrium\n');
disp('A = ');
disp(A)
disp('B = ');
disp(B)
disp('C = ');
disp(Cout)
disp('D = ');
disp(D)

%% PART 3: Open-loop transfer function
sys_ss = ss(A, B, Cout, D);
H = tf(sys_ss);

fprintf('\nPart 3: Open-Loop Transfer Function H(s) = delta_vdc / delta_vs\n');
disp(H)

%% PART 3 & 5: Stability check
% Stability condition from project:
%   P / vdc0^2 < 1 / R
lhs = P / vdc0^2;
rhs = 1 / R;

traceA = trace(A);
detA   = det(A);
eigA   = eig(A);

fprintf('\nPart 5: Stability Check\n');
fprintf('P/vdc0^2 = %.6f\n', lhs);
fprintf('1/R      = %.6f\n', rhs);
fprintf('trace(A) = %.6f\n', traceA);
fprintf('det(A)   = %.6f\n', detA);
disp('eig(A) = ');
disp(eigA)

if lhs < rhs && all(real(eigA) < 0)
    fprintf('Result: Linearized system is stable.\n');
else
    fprintf('Result: Linearized system is NOT stable.\n');
end

%% PART 6: Controller design
%   Use PD controller
%   Kp = 0
%   Kd = 2.5e-4
%
Kp = 0;
Ki = 0;
Kd = 2.5e-4;

s = tf('s');
Gc = Kp + Kd*s;   

T_cl = feedback(Gc * H, 1); 

fprintf('\nPart 6: PD Controller\n');
fprintf('Kp = %.6g\n', Kp);
fprintf('Ki = %.6g\n', Ki);
fprintf('Kd = %.6g\n', Kd);

fprintf('\nClosed-loop transfer function:\n');
disp(T_cl)

cl_poles = pole(T_cl);
fprintf('Closed-loop poles:\n');
disp(cl_poles)

step_info = stepinfo(T_cl);
fprintf('Step response info:\n');
disp(step_info)

%% Plot linear closed-loop step response
figure;
step(T_cl, 0:1e-5:0.01);
grid on;
title('Linear Closed-Loop Step Response with PD Controller');
xlabel('Time (s)');
ylabel('\delta v_{dc}');

%% PART 6: Linearized response from initial condition delta x(0) = [8; -8]
dx0 = [8; -8];

t_lin = linspace(0, 0.01, 2000);
[y_lin_init, t_lin_init, x_lin_init] = initial(sys_ss, dx0, t_lin);

figure;
plot(t_lin_init, y_lin_init, 'LineWidth', 1.5);
grid on;
title('Linearized Output Response from Initial Condition (Open-Loop State-Space)');
xlabel('Time (s)');
ylabel('\delta v_{dc}');

%% PART 7: Nonlinear closed-loop simulation with safety bounds
x0_nl = [iL0; vdc0] + dx0;   % [8.33; 72]

fprintf('\nPart 7: Nonlinear Closed-Loop Simulation\n');
fprintf('Initial nonlinear state:\n');
fprintf('iL(0)  = %.4f A\n', x0_nl(1));
fprintf('vdc(0) = %.4f V\n', x0_nl(2));

% Simulation time
tspan = [0 0.01];

opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

[t_nl, x_nl] = ode23s(@(t,x) nonlinear_closed_loop_ode(t, x, R, P, L, C, vs0, Kp, Kd), ...
                      tspan, x0_nl, opts);

iL_nl  = x_nl(:,1);
vdc_nl = x_nl(:,2);

% Reconstruct controller input and actual source voltage
vs_nl = zeros(size(t_nl));
delta_vs_nl = zeros(size(t_nl));
for k = 1:length(t_nl)
    xk = x_nl(k,:)';
    [~, vs_temp, delta_vs_temp] = nonlinear_closed_loop_ode(t_nl(k), xk, R, P, L, C, vs0, Kp, Kd);
    vs_nl(k) = vs_temp;
    delta_vs_nl(k) = delta_vs_temp;
end

%% Nonlinear plots
figure;
plot(t_nl, vdc_nl, 'LineWidth', 1.6);
hold on;
yline(80, '--');
grid on;
title('Nonlinear Closed-Loop Capacitor Voltage');
xlabel('Time (s)');
ylabel('v_{dc} (V)');
legend('v_{dc}(t)','80 V reference','Location','best');

figure;
plot(t_nl, iL_nl, 'LineWidth', 1.6);
grid on;
title('Nonlinear Closed-Loop Inductor Current');
xlabel('Time (s)');
ylabel('i_L (A)');

figure;
plot(t_nl, vs_nl, 'LineWidth', 1.6);
grid on;
title('Nonlinear Closed-Loop Source Voltage');
xlabel('Time (s)');
ylabel('v_s (V)');

%% Basic nonlinear performance summary
vdc_error = abs(vdc_nl - 80);
tol = 0.02 * 80;  % 2% band around 80 V
idx_settle = find(vdc_error <= tol, 1, 'first');

if isempty(idx_settle)
    Ts_nl = NaN;
else
    % Make sure it stays there afterward
    Ts_nl = NaN;
    for k = idx_settle:length(t_nl)
        if all(vdc_error(k:end) <= tol)
            Ts_nl = t_nl(k);
            break;
        end
    end
end

overshoot_nl = max(vdc_nl - 80);

fprintf('\nNonlinear Closed-Loop Results\n');
fprintf('Max iL    = %.4f A\n', max(iL_nl));
fprintf('Max vdc   = %.4f V\n', max(vdc_nl));
fprintf('Max vs    = %.4f V\n', max(vs_nl));
fprintf('Overshoot above 80 V = %.6f V\n', overshoot_nl);
fprintf('Approx. 2%% settling time = %.6e s\n', Ts_nl);

if all(vs_nl >= 0 & vs_nl <= 160) && all(iL_nl >= 0 & iL_nl <= 20) && all(vdc_nl >= 0 & vdc_nl <= 200)
    fprintf('Safety bounds satisfied during simulation.\n');
else
    fprintf('Warning: One or more safety bounds were violated.\n');
end

fprintf('\nDone.\n');

%% Local function: nonlinear closed-loop system
function [xdot, vs, delta_vs] = nonlinear_closed_loop_ode(~, x, R, P, L, C, vs0, Kp, Kd)

    % Unpack states
    iL  = x(1);
    vdc = x(2);

    % Safety clamp on states
    iL  = min(max(iL, 0), 20);
    vdc = min(max(vdc, 1e-3), 200);   % small positive lower bound to avoid division by zero

    % Error for nonlinear implementation
    e = vs0 - vdc;   % target is 80 V

    % Nonlinear capacitor dynamics used to estimate de/dt
    dvdc_dt_open = (1/C)*(iL - vdc/R - P/vdc);

    % Since e = vs0 - vdc, then de/dt = -dvdc/dt
    de_dt = -dvdc_dt_open;

    % PD controller output in deviation variables
    delta_vs = Kp*e + Kd*de_dt;

    % Translate back to actual plant input
    vs = vs0 + delta_vs;

    % Source voltage safety saturation
    vs = min(max(vs, 0), 160);

    % Nonlinear plant dynamics
    diL_dt  = (1/L)*(vs - vdc);
    dvdc_dt = (1/C)*(iL - vdc/R - P/vdc);

    
    if (x(1) <= 0 && diL_dt < 0) || (x(1) >= 20 && diL_dt > 0)
        diL_dt = 0;
    end
    if (x(2) <= 0 && dvdc_dt < 0) || (x(2) >= 200 && dvdc_dt > 0)
        dvdc_dt = 0;
    end

    xdot = [diL_dt; dvdc_dt];
end