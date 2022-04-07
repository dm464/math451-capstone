function UU = funnel_r_time_evolution2(u_pre)
%% Time evolution of a thin film of liquid on a funnel - constant volume
% Domain: r \in [R, L]
% Input: precursor film 'u_pre'.
% NOTE!! 
% The setup should be EXACTLY the same as those in 'funnel_r_precursor.m'
% Te-Sheng Lin
%% Parameters
% funnel angle
alpha = 47*pi/180;
% number of grid points
N = 500;
G=1;
bb = 0.1;
%%
SA = sin(alpha);
CA = cos(alpha);
D = 1;
%% domain r\in[R, L]
R = 1;
L = 100;
%% Grid construction
% spatial grid step size
dx = (L-R)/N;
% grid points
x = R + dx/2 + (0:N-1)'*dx;
%% initial guess/condition
u0 = u_pre +((1 - tanh(2*(x-(L-10)))).*(1 + tanh(2*(x-(L-30))))*(1-bb)/4);
% u0 = u_pre +((1 - tanh(2*(r-(L-40)))).*(1 + tanh(2*(r-(L-60))))*(1-bb)/4);
% sub = 0.1*ones(size(r));
sub = 0.2*((1 - tanh(2*(x-(L-60)))).*(1 + tanh(2*(x-(L-70))))*(1-bb)/4);

% u0 = u_pre + 3.83.*((1 - tanh(2*(r-(L-10)))).*(1 + tanh(2*(r-(L-40))))*(1-bb)/4);
%%
dt = 1;
it_num = 100;
time = zeros(it_num+1, 1);
measure = zeros(it_num+1, 4);
measure(1,:) = [80, bb, 80, 1];
UU=zeros(it_num+1,N); %%
UU(1,:)=u0; %%
%%
%figure
%titlehandle = suptitle(['time=', num2str(0, '%7.1f')]);
% plot([R*CA, L*CA], [R*SA, L*SA], 'b', 'LineWidth', 2)
% hold on
% h1 = plot(r*CA-SA*u0, r*SA+CA*u0, 'r', 'LineWidth', 2);
% xlabel('x'); ylabel('z');
%drawnow
%%
figure(1)
%titlehandle = suptitle(['time=', num2str(0, '%7.1f')]);
plot(x, sub, x, u0+sub');
drawnow
%%
continue_flag = 1;
ii = 1;
while continue_flag
    time(ii+1) = ii*dt;
    [~, y_all]=ode15s(@myfun,[ii-1 ii]*dt, u0);
    y = y_all(end,:)';
    ii*dt;
    %h1.XData = r*CA-SA*y;
    %h1.YData = r*SA+CA*y;
    %set(titlehandle, 'string', ['time=', num2str(time(ii+1), '%7.1f')]);
    drawnow
    %% Some measure
    % measure(:,1:2): location and height of the min.
    [~, I2] = min(y-u_pre);
    measure(ii+1,1:2) = [x(I2), y(I2)];
    % measure(:,3:4): location and height of the max.
    [~, I2] = max(y-u_pre);
    measure(ii+1,3:4) = [x(I2), y(I2)];
    %if ii>10
    %    continue_flag=0;
    %end
    if (measure(ii,3)<10) && ii>10
        continue_flag=0;
    end
    u0 = y_all(end, :);
    fluid_on_sub = u0+sub';
    plot(x, sub, x, fluid_on_sub);
    legend('Substrate', 'Fluid');
    ylim([0 1.5])
    drawnow;
    ii = ii+1;
    UU(ii,:)=u0; %%
    %ii
end
%%
figure(2)
time = time(1:ii, :);
measure = measure(1:ii, :);
subplot(3,1,1)
plot(time, measure(:,1), time, measure(:,3), 'LineWidth', 2)
title('Time evolution of locations at the max. and min.');
xlabel('time');
ylabel('r');
subplot(3,1,2)
plot(time, measure(:,2), 'LineWidth', 2)
title('Time evolution of the min. height');
xlabel('time');
subplot(3,1,3)
plot(time, measure(:,4), 'LineWidth', 2)
title('Time evolution of the max. height');
xlabel('time');

function F = myfun(~,h)
    sh = h+sub;
    %% Ghost points outside the computational domain
    % hn1 = h(R-dr/2)
    % hn2 = h(R-3*dr/2)
    % hp1 = h(L+dr/2)
    % hp2 = h(L+3*dr/2)
    % NOTE!!
    % The boundary conditions related to hn2 and hp2 are imposed later, 
    % so these two points need not to be defined here, 
    % we simply assign 0 to them. 
    h(1) = u_pre; h(N) = u_pre;
    hn1 = h(1);
    hn2 = u_pre;
    %hp1 = 2*bb-h(N);
    hp1 = u_pre;
    hp2 = u_pre;
    shl = sh(1);
    shr = sh(N);
    %% h3: Evaluate h^3 at N+1 regular grid points
    h3 = ([h.^3;hp1^3]+[hn1^3;h.^3])/2;
    
    %% tmp1: Evaluate (s+h)_x_x at N+2 grid points (central difference)
    FDpos = ([sh;shr;shr]-[shl;sh;shr])/dx; %Forward diff
    FDpos(2)=0;
    FDpos(N+1)=0; % TODO change boundary conditions
    FDneg = ([shl;sh;shr]-[shl;shl;sh])/dx; %Forward diff
    FDneg(3)=0;
    FDneg(N+2)=0;
    tmp1 = (FDpos-FDneg)/dx; %Backward diff
    % tmp1 = ([h;hp1;hp2] - 2*[hn1;h;hp1] + [hn2;hn1;h])/(dr^2);
    
    %% tmp2: Evaluate ((1/r)(r*h_r)_r)_r -TA/r^2 - G*(SA+CA*h_r) at N+1 grid points    
    % NOTE!!
    % Some boundary conditions are imposed here
    % Evaluate h_x_x_x at N+1 grid points (forward difference 0-N)
    tmp2 = (tmp1(2:(N+2)) - tmp1(1:(N+1)))/dx;
    % Boundary condition: ((1/r)(r*h_r)_r)_r=0 at r=R
    tmp2(1)=0;
    % Evaluate h_x_x_x - D*h_x at N+1 grid points (forward difference)
    tmp2 = tmp2 - D*CA*([sh;shr]-[shl;sh])/dx;
    % Boundary condition: ((1/r)(r*h_r)_r)_r - G*CA*h_r=0 at r=L
    tmp2(N+1)=0;    
    % Evaluate h_x_x_x - D*h_x - 1 at N+1 grid points
    tmp2 = tmp2 - SA;
    %% Q: Evaluate h3*(h_x_x_x - D*h_x - 1) at N+1 grid points
    % q = h3*tmp2 is the flux
    Q = h3.*tmp2;
    %% F = Q_x (backward difference)
    F = -(Q(2:N+1)-Q(1:N))./(dx);
    % F(9*N/10)=0.0;
end
end