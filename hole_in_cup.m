funnel_r_time_evolution(1e-1);
function UU = funnel_r_time_evolution(u_pre)
%% Time evolution of a thin film of liquid on a funnel - constant volume
%% Parameters
alpha = 0;
% number of grid points
N = 500;
bb = 0.1;
%%
SA = sin(alpha);
CA = cos(alpha);
hc = 200e-6;
rho = 900;
g = 9.8;
gamma = 20e-3;
mu = 50e-6;
xc = sqrt(gamma/(rho*g));
tc = (3*mu*xc)/gamma;
Ca = (hc^2*rho*g)/(3*gamma);
D = (3*Ca)^(1/3);
%% domain r\in[R, L]
R = 0;
L = 100;
%% Grid construction
% spatial grid step size
dx = (L-R)/N;
% grid points
x = R + dx/2 + (0:N-1)'*dx;
%% initial guess/condition
u0 = 0.9+u_pre -((1 - tanh(2*(x-(L-30)))).*(1 + tanh(2*(x-(L-70))))*(1-bb)/4);
true_area = 64;
u_init = u0;
sub = zeros(size(u0));
u02=u0;

%%
epsilons = 1.5;%0.5:0.25:2;
time_to_close = zeros(size(epsilons));
time_to_touch = zeros(size(epsilons));
tol = 0.05; % 0.01
dt = 0.1; %0.1
it_num = 100;
time = zeros(it_num+1, 1);
%%
figure(1)
plot(x,u02);
ylim([0 1]);
legend('Fluid');
drawnow;
%%
for j=1:length(epsilons)
    u0 = u_init;
    u02 = u_init;
    continue_flag = 1;
    ii = 1;
    while continue_flag
        current_area = trapz(x, u02);
        time(ii+1) = ii*dt;
        eps=epsilons(j);
        om=10;
        G = nthroot(1+eps*cos(om*time(ii+1)),3);
        V = G*D*CA;
        [~, y_all]=ode15s(@myfun,[ii-1 ii]*dt, u0);
        y = y_all(end,:)';
        drawnow;
        %% Some measure
        [~, mIndex] = min(u02);
        if time_to_touch(j) == 0 && (mIndex == (N/2) || mIndex == (N/2+1))
            time_to_touch(j) = time(ii+1);
        end
        if max(u02)-min(u02) <= tol
            continue_flag=0;
            time_to_close(j) = time(ii+1);
        end
        if abs(current_area - true_area) > (tol * true_area)
            continue_flag=0;
            fprintf('The area is no longer consistent');
        end
        u0 = y_all(end, :);
        u02=u0+sub';
        plot(x,u02);
        ylim([0 1.2])
        title(['\epsilon = ', num2str(eps)]);
        legend('Fluid');
        drawnow;
        ii = ii+1;
    end
end
time_to_close
time_to_touch
scatter(epsilons, time_to_touch)
hold on
scatter(epsilons, time_to_close)
hold off
title(['\epsilon vs time (tol = ', num2str(tol), ')'])
xlabel('\epsilon')
ylabel('time (t_s)')
legend('Time to merge', 'Time to close')

function F = myfun(~,h)
    %% Ghost points outside the computational domain
    % hn1 = h(R-dr/2)
    % hn2 = h(R-3*dr/2)
    % hp1 = h(L+dr/2)
    % hp2 = h(L+3*dr/2)
    % NOTE!!
    % The boundary conditions related to hn2 and hp2 are imposed later, 
    % so these two points need not to be defined here, 
    % we simply assign 0 to them. 
    sh=h+sub;
    shl1 = sh(1);
    shl2 = sh(2);
    shr1 = sh(N);
    shr2 = sh(N-1);
    hn1 = h(1);
    hp1 = h(N);
    %% h3: Evaluate h^3 at N+1 regular grid points
    h3 = ([h.^3;hp1^3]+[hn1^3;h.^3])/2;
    %% tmp1: Evaluate h_x_x at N+2 grid points (central difference)
    % tmp1 = ([h;hp1;hp2] - 2*[hn1;h;hp1] + [hn2;hn1;h])/(dr^2);  %central
    FDpos=([sh;shr1;shr2]-[shl1;sh;shr1])/dx; %backward difference with (1 to 502)
    FDneg=([shl1;sh;shr1]-[shl2;shl1;sh])/dx; %backward difference with (0 to 501)
    tmp1=(FDpos-FDneg)/dx; %forward difference with (0 to 501) (same as central)
    %% tmp2: Evaluate ((1/r)(r*h_r)_r)_r -TA/r^2 - G*(SA+CA*h_r) at N+1 grid points    
    % NOTE!!
    % Some boundary conditions are imposed here
    % Evaluate h_x_x_x at N+1 grid points (backward difference 1-501)
    tmp2 = (tmp1(2:(N+2)) - tmp1(1:(N+1)))/dx;
    % Boundary condition: h_x_x_x=0 at r=R=L
    tmp2(1)=0;
    tmp2(N)=0;
    % Evaluate h_x_x_x - D*h_x at N+1 grid points
    h_x=FDpos(1:501);
    %tmp2 = tmp2 - G*D*CA*([sh;shr1]-[shl1;sh])/dx;
    tmp2 = tmp2 - V*(h_x)/dx; %3rd derivative
    % Boundary condition: ((1/r)(r*h_r)_r)_r - G*CA*h_r=0 at r=L
    % Evaluate h_x_x_x - D*h_x - 1 at N+1 grid points
    tmp2 = tmp2 - SA;
    %% Q: Evaluate h3*(h_x_x_x - D*h_x - 1) at N+1 grid points
    % q = h3*tmp2 is the flux
    Q = h3.*tmp2; 
    %% F = Q_x (forward difference)
    F = -(Q(2:N+1)-Q(1:N))./(dx);
    %F(9*N/10)=0;
end
end