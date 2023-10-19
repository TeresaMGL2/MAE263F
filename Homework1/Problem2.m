clear all
NN = [3,11,21,51,101];% Number of nodes
dtt = [10,1,1e-1,5e-2,1e-2]; % second - Time step size
final_v_n = zeros(5);
final_v_t = zeros(5);

%-----------------------change number of nodes-------------
% time step unchanged = 1e-2
for cc = 1:length(dtt)
N = NN(cc);
dt = 1e-2;
ndof = N * 2; % number of degrees of freedom

L = 0.1; %
dL = L / (N-1);

% Radii of spheres
R = zeros(N,1); % Vector of size N - Radius of N nodes
midNode = (N+1)/2;
R(:) = dL/10;
R(midNode) = 0.025;

% Density
rho_m = 7000; % kg/m^3
rho_f = 1000; % fluid
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
T = 50; % second - total simulation time

% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - edge bending stiffness
EA = Y * pi * r0^2; % Newton

% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * dL; % x coordinates
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_m; % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Viscous damping matrix, C
C = zeros(ndof,ndof);
for k=1:N
    C(2*k-1, 2*k-1) = 6 * pi * visc * R(k); %viscous damping about x
    C(2*k, 2*k) = C(2*k-1, 2*k-1); %viscous damping about y
end

% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
    W(2*k-1) = 0; % weight along x is zero
    W(2*k) = -4/3*pi*R(k)^3*(rho_m-rho_f)*g;
end

% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end
u0 = zeros(ndof, 1); % old velocity (initial velocity)
% tolerance
tol = EI/L^2 * 1e-3; % small enouch force that can be neglected
% Time marching scheme
Nsteps = round(T/dt);
% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);

for c = 2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt);
    fprintf('number of nodes = %f\n', N);
    % Guess
    q = q0; % New DOFs are initialized to be equal to old DOFs
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        f = M / dt * ( (q-q0)/dt - u0 );
        J = M / dt^2;
        %
        % Elastic forces
        % -Linear spring
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            l_k = dL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end
        % -Bending spring
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            l_k = dL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end

        % Viscous force
        f = f + C * (q-q0) / dt;
        J = J + C / dt;

        % Weight
        f = f - W;

        % At this point, we have f and J
        % Update
        dq = J \ f;
        q = q - dq;
        err = sum ( abs(f) );
    end
    % New velocity
    u = (q - q0) / dt;
    % Store some information
    all_mid_v(c) = u(2*midNode); %stores y-velo of central node
    all_mid_x(c) = q(2*midNode); %stores y-velo of central node

%     % Plot
%     figure(1);
%     plot( q(1:2:end), q(2:2:end), 'ro-');
%     axis equal
%     xlabel('x [meter]');
%     ylabel('y [meter]');
%     drawnow
    % Update (new becomes old)
    q0 = q;
    u0 = u;
end

final_v_n(cc) = u(end);
% figure(1);
% plot( q(1:2:end), q(2:2:end), 'ro-');
% axis equal
% xlabel('x [meter]');
% ylabel('y [meter]');



end

% Plot middle node downward velocity
figure(1);
plot(NN, final_v_n, 'k-');
xlabel('Number of Nodes');
ylabel('Number of Nodes vs. Terminal Velocity [m/s]');


%---------------------change time step-----------------
% number of nodes unchanged = 21
for cc = 1:length(dtt)
N = 21;
dt = dtt(cc);
ndof = N * 2; % number of degrees of freedom

L = 0.1; %
dL = L / (N-1);

% Radii of spheres
R = zeros(N,1); % Vector of size N - Radius of N nodes
midNode = (N+1)/2;
R(:) = dL/10;
R(midNode) = 0.025;

% Density
rho_m = 7000; % kg/m^3
rho_f = 1000; % fluid
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
T = 50; % second - total simulation time

% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - edge bending stiffness
EA = Y * pi * r0^2; % Newton

% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * dL; % x coordinates
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_m; % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Viscous damping matrix, C
C = zeros(ndof,ndof);
for k=1:N
    C(2*k-1, 2*k-1) = 6 * pi * visc * R(k); %viscous damping about x
    C(2*k, 2*k) = C(2*k-1, 2*k-1); %viscous damping about y
end

% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
    W(2*k-1) = 0; % weight along x is zero
    W(2*k) = -4/3*pi*R(k)^3*(rho_m-rho_f)*g;
end

% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end
u0 = zeros(ndof, 1); % old velocity (initial velocity)
% tolerance
tol = EI/L^2 * 1e-3; % small enouch force that can be neglected
% Time marching scheme
Nsteps = round(T/dt);
% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);

for c = 2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt);
    fprintf('Time step = %f\n', dt);
    % Guess
    q = q0; % New DOFs are initialized to be equal to old DOFs
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        f = M / dt * ( (q-q0)/dt - u0 );
        J = M / dt^2;
        %
        % Elastic forces
        % -Linear spring
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            l_k = dL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end
        % -Bending spring
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            l_k = dL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end

        % Viscous force
        f = f + C * (q-q0) / dt;
        J = J + C / dt;

        % Weight
        f = f - W;

        % At this point, we have f and J
        % Update
        dq = J \ f;
        q = q - dq;
        err = sum ( abs(f) );
    end
    % New velocity
    u = (q - q0) / dt;
    % Store some information
    all_mid_v(c) = u(2*midNode); %stores y-velo of central node
    all_mid_x(c) = q(2*midNode); %stores y-velo of central node

%     % Plot
%     figure(1);
%     plot( q(1:2:end), q(2:2:end), 'ro-');
%     axis equal
%     xlabel('x [meter]');
%     ylabel('y [meter]');
%     drawnow
    % Update (new becomes old)
    q0 = q;
    u0 = u;
end

final_v_t(cc) = u(end);
% figure(1);
% plot( q(1:2:end), q(2:2:end), 'ro-');
% axis equal
% xlabel('x [meter]');
% ylabel('y [meter]');



end

% Plot middle node downward velocity
figure(2);
plot(dtt, final_v_t, 'k-');
xlabel('Time Step [s]');
ylabel('Time Step [s] vs. Terminal Velocity [m/s]');

%----------------------borrowed functions
function F = gradEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the derivative of stretching energy E_k^s with 
% respect to x_{k-1}, y_{k-1}, x_k, and y_k.
%
F = zeros(4,1);

F(1) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk);
F(2) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk);
F(3) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * xkp1 - 0.2e1 * xk);
F(4) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * ykp1 - 0.2e1 * yk);

F = 0.5 * EA * l_k * F;

end

function J = hessEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the 4x4 hessian of the stretching energy E_k^s with
% respect to x_k, y_k, x_{k+1}, and y_{k+1}.
%
J11 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * xkp1 + 2 * xk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((-2 * xkp1 + 2 * xk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J12 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * ykp1 + 2 * yk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (-2 * ykp1 + 2 * yk) / 0.2e1;
J13 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (2 * xkp1 - 2 * xk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J14 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (2 * ykp1 - 2 * yk) / 0.2e1;
J22 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * ykp1 + 2 * yk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((-2 * ykp1 + 2 * yk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J23 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) * (-2 * ykp1 + 2 * yk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * ykp1 + 2 * yk) * (2 * xkp1 - 2 * xk) / 0.2e1;
J24 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (-2 * ykp1 + 2 * yk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * ykp1 + 2 * yk) * (2 * ykp1 - 2 * yk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J33 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((2 * xkp1 - 2 * xk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J34 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (2 * xkp1 - 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (2 * xkp1 - 2 * xk) * (2 * ykp1 - 2 * yk) / 0.2e1;
J44 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((2 * ykp1 - 2 * yk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;

J = [J11 J12 J13 J14;
     J12 J22 J23 J24;
     J13 J23 J33 J34;
     J14 J24 J34 J44];

J = 0.5 * EA * l_k * J;

end

function dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI)

%
% This function returns the derivative of bending energy E_k^b with respect
% to x_{k-1}, y_{k-1}, x_k, y_k, x_{k+1}, and y_{k+1}.
%
% curvature0 is the "discrete" natural curvature [dimensionless] at node (xk, yk).
% l_k is the voronoi length of node (xk, yk).
% EI is the bending stiffness.
%

node0 = [xkm1, ykm1, 0];
node1 = [xk, yk, 0];
node2 = [xkp1, ykp1, 0];
%     m1e, 
m2e = [0 0 1];
%     m1f,
m2f = [0 0 1];
kappaBar = curvature0;

%% Computation of gradient of the two curvatures
gradKappa = zeros(6,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d2 = (m2e + m2f) / chi;

% Curvatures
kappa1 = kb(3); % 0.5 * dot( kb, m2e + m2f); % CHECKED

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

gradKappa(1:2, 1) = -Dkappa1De(1:2);
gradKappa(3:4, 1) = Dkappa1De(1:2) - Dkappa1Df(1:2);
gradKappa(5:6, 1) = Dkappa1Df(1:2);

%% Gradient of Eb
dkappa = kappa1 - kappaBar;
dF = gradKappa * EI * dkappa / l_k;
end

