%% Parameters
% funnel angle
u_pre=1e-2;
alpha = 47*pi/180;
% number of grid points
N = 500;
% 
G=1;
 
bb = 0.1;
 
%%
SA = sin(alpha);
CA = cos(alpha);
TA = tan(alpha);
 
%% domain r\in[R, L]
R = 1;
L = 100;
 
%% Grid construction
% spatial grid step size
dr = (L-R)/N;
% grid points
r = R + dr/2 + (0:N-1)'*dr;
 
% r2: N+2 points including ghost points
r2 = [r(1)-dr;r;r(N)+dr];
% r3: N+1 points on regular grid
r3 = R + (0:N)'*dr;
 
%% initial guess
u0 = u_pre +((1 - tanh(2*(r-(L-10)))).*(1 + tanh(2*(r-(L-30))))*(1-bb)/4);
%u0 = u_pre + 3.83.*((1 - tanh(2*(r-(L-10)))).*(1 + tanh(2*(r-(L-40))))*(1-bb)/4);

plot(r,u0)