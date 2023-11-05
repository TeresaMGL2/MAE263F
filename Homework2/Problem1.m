clear all


global Fg M dt
global kappaBar EI GJ voronoiLength
global EA refL

%% Parameters
%Physical parameters
l = 0.2; %rod length
rn = 0.02; %natural curve
rho=1000;
r0=0.001; %rod radius
Y = 10e6; %Young's
G = Y/3; %shear modulus
g = -9.81; %gravitational constant

% Stifness variables
EI = Y * pi * r0^4 / 4; % Bending stiffness
GJ = G * pi * r0^4/2; % Shearing stiffness
EA = Y * pi * r0^2; % Stretching stiffness

%Partitions
n = 50; %number of nodes
dt = 0.01; %time step
dtheta = (l/rn)*(1/(n-1));
T=5; %total time of simulation
fixed = 1:7; %index for two nodes and first angle
free = 8:4*n-1;
tol=1/Y;

%% Initializing q
q = zeros(4*n-1,1); %stores [x1,y1,z1, theta1 to xn,yn,zn]
k=1;
for i=1:4:4*n %compute initial rod position
    q(i) = rn*cos((k-1)*dtheta);
    q(i+1) = rn*sin((k-1)*dtheta);
    k=k+1;
end
q0 = q;
u = zeros(4*n-1,1);

%% get reference length 
refL = zeros(n-1, 1);
for i=1:n-1 % loop over each edge
 dx = q(i*4+1:i*4+3) - q(i*4-3:i*4-1); %distance vector from xi to xi-1
 refL(i) = norm(dx);
end

%% Voronoi length (length associated with each node; used for bending and twisting)
voronoiLength = zeros(n, 1);
for c=1:n
 if c==1
 voronoiLength(c) = 1/2 * refL(c);
 elseif c==n
 voronoiLength(c) = 1/2 * refL(c-1);
 else
 voronoiLength(c) = 1/2 * refL(c-1) + 1/2 * refL(c);
 end
end

%% initialize time parallel frame 
a1 = zeros(n-1, 3); % First reference director
a2 = zeros(n-1, 3); 
tangent = computeTangent(q); % calculate tangent for all edges

t0 = tangent(1, :); % tangent on the first edge
t1 = [0;0;1]; % arbitrary start
a1Tmp = cross(t0, t1); 
if abs(a1Tmp) < 1e-6 % if overlay too large
 t1 = [1;0;0]; %change the arbitrary start
 a1Tmp = cross(t0, t1);
end
a1(1,:) = a1Tmp / norm(a1Tmp); % create the initial a1
a2(1,:) = cross(t0, a1(1,:)); %get a2 from a1 and t0
% parallel transport to get the rest
for i=2:n-1
 t0 = tangent(i-1,:); % previous tangent on c-1-th edge
 t1 = tangent(i,:); % current tangent on c-th edge
 a1_0 = a1(i-1,:); % previous a1 director on c-1-th edge
 a1_1 = parallel_transport( a1_0, t0, t1);
 a1(i,:) = a1_1 / norm( a1_1 );
 a2(i,:) = cross(t1, a1(i,:));
end

%% get material frame
theta = q(4:4:4*n-3); 
[m1, m2] = computeMaterialDirectors(a1, a2, theta);

%% Reference twist
refTwist = zeros(n, 1);

%% Natural curvature
kappaBar = getkappa( q, m1, m2 ); % Natural curvature at each node

%% weight assignment
% given rod, only edge has weight
Mv = zeros(4*n-1, 1);
totalM = pi*r0^2*l*rho; 
dm = totalM / (n-1); % mass per edge
for i=1:n
 if i==1 || i==n
 Mv(i*4-3) = dm/2;
 Mv(i*4-2) = dm/2;
 Mv(i*4-1) = dm/2;
 else
 Mv(i*4-3) = dm;
 Mv(i*4-2) = dm;
 Mv(i*4-1) = dm;
 end
end
for c=1:n-1 % Loop over edges
Mv(4*c) = 1/2 * dm * r0^2;
end
M = diag(Mv); % ndof x ndof sized mass matrix

Fg = zeros(4*n-1, 1);
for i=3:4:4*n-1
    Fg(i) =Mv(i)*g;
end

%% Numerical steps
for k = 1:T/dt
 fprintf('Current time = %f\n', k*dt);
 [q, u, a1, a2] = objfun(q0, u, a1, a2, free, tol, refTwist);


 % store previous position
 q0 = q;
 

 % Store z of last node
 z(k) = q(length(q));

 %plot initial orientation of rods
 if k==1
    figure(1)
    theta = q(4:4:end);
    [m1, m2] = computeMaterialDirectors(a1, a2, theta);
    plotrod(q, a1, a2, m1, m2, 0);
 end

 %Plot
 if mod(k,200) == 0 %plot every 200 time step
    theta = q(4:4:end);
    [m1, m2] = computeMaterialDirectors(a1, a2, theta);
    plotrod(q, a1, a2, m1, m2, k*dt);
 end

 %plot end orientation of rods
 if k==T/dt
    theta = q(4:4:end);
    [m1, m2] = computeMaterialDirectors(a1, a2, theta);
    plotrod(q, a1, a2, m1, m2, T);
 end
end

figure(2)
%% Plotting the z position of the last node across time
plot(0:dt:T-dt,z,'*--') 
title('Z position of last node across time')
xlabel('Time(s)')
ylabel('Z position(m)')
ylim([-0.08,0])
xlim([0,5])

function [q, u, a1, a2] = objfun(q0, u, a1, a2, free, ...
 tol, refTwist)
global M dt Fg

% initial guess
q = q0;
iter = 1;
error = 1;
while error > tol

 % Compute reference frame
 [a1it, a2it] = computeTimeParallel(a1, q0, q);

 % Compute reference twist
 tangent = computeTangent(q);
 refTwist_it = computeRefTwist(a1it, tangent, refTwist);

 % Compute material frame
 theta = q(4:4:length(q));
 [m1, m2] = computeMaterialDirectors(a1it, a2it, theta);

 % Compute elastic force and Jacobian
 [Fb, Jb] = getFb(q, m1, m2); % Bending
 [Ft, Jt] = getFt(q, refTwist_it); % Twisting
 [Fs, Js] = getFs(q); % Stretching
 % Equations of motion
 f = M/dt*((q-q0)/dt - u) -Fb -Ft -Fs -Fg;
 J = M/dt^2 - Jb - Jt - Js;

 f_free = f(free);
 J_free = J(free, free);

 % Newton's update
 dq_free = J_free \ f_free ;
 q(free) = q(free) - dq_free;

 % Error
 error = sum( abs( f_free ) );

 fprintf('Iter = %d, error=%f\n', iter, error);
 iter = iter + 1;

end
u = (q - q0) / dt;
a1 = a1it;

a2 = a2it;
end



function [m1,m2] = computeMaterialDirectors(a1, a2, theta)
% Inputs:
% a1 is a matrix of size ne x 3. This contains the first
% reference director on each edge.
% a2 is a matrix of size ne x 3.
% theta is a vector of size ne ( theta = q(4:4:end) )
%
% Outputs:
% m1 and m2 are matrices of size ne x 3. Each column of
% these matrices contain the material directors on each
% edge.
ne = length(theta);
m1 = zeros(ne,3);
m2 = zeros(ne,3);
for c=1:ne % loop over edges
 cs = cos(theta(c));
 ss = sin(theta(c));
 m1(c,:) = cs * a1(c,:) + ss * a2(c,:);
 m2(c,:) = - ss * a1(c,:) + cs * a2(c,:);
end
end

function refTwist = computeRefTwist(a1, tangent, refTwist)
% a1 is a matrix of size ne x 3. Each column of a1
% contains the first time-parallel reference director
% on each edge.
%
% tangent is a matrix of size ne x 3. It contains all
% the tangent vectors on the edges.
%
% refTwist is a vector of size nv (one reference twist
% for each node). This (optionally) is provided as an
% input (guess solution).
[ne, ~] = size(a1); % number of edges
nv = ne + 1;
% refTwist = zeros(nv, 1);
for c=2:ne % all internal nodes

 u0 = a1(c-1,:); % a1 vector of previous edge
 u1 = a1(c, :); % a1 vector of current edge
 t0 = tangent(c-1,:); % tangent of previous edge
 t1 = tangent(c,:); % tangent of current edge
 ut = parallel_transport(u0, t0, t1);
 % ut and u1 are the same?

 % Method 1: Okay? But 2*pi issue?
 % refTwist(c) = signedAngle(ut, u1, t1);

 % Method 2
 ut = rotateAxisAngle( ut, t1, refTwist(c) );
 refTwist(c) = refTwist(c) + signedAngle(ut, u1, t1);
end
end

function tangent = computeTangent(q)
% tangent is a matrix of size ne x 3
% q is a vector of size 4*nv-1
nv = (length(q)+1)/4;
ne = nv - 1;
tangent = zeros(ne, 3);
% Each column of "tangent" matrix is a 3-element
% vector representing the tangent on each edge
for c = 1:ne
 xc = q(4*c-3: 4*c-1);
 xcp1 = q(4*c+1: 4*c+3);
 edge = xcp1 - xc;
 tangent(c,:) = edge / norm(edge);
end
end

function [a1, a2] = computeTimeParallel(a1_old, q0, q)
% a1_old is a matrix of size nex3: It contains the
% time-parallel reference director (a1) at each of
% the ne edges at the old time step
% q0 is the DOF vector at the old time step
% q is the DOF vector at the new time step
% Outputs:
% a1 is the time-parallel reference director (a1) at
% each of the ne edges at the new time step. Size of
% a1 is ne x 3.
% a2 is the second reference director. It's size is
% ne x 3.
nv = (length(q)+1)/4;
ne = nv - 1;
tangent0 = computeTangent(q0); % Tangents at old step
tangent = computeTangent(q); % Tangents at new step
a1 = zeros(ne, 3);
a2 = zeros(ne, 3);
for c=1:ne % loop over edges
 t0 = tangent0(c,:); % tangent on c-th edge at old step
 t = tangent(c,:); % tangent on c-th edge at new step

 a1_local_old = a1_old(c,:); % a1 director on c-th edge
 % at old step
 a1_local = parallel_transport( a1_local_old, t0, t);
 % a1_local is the first reference director on c-th
 % edge at new step

 % Just to be careful: enforce a1 and t are perp.
 a1_local = a1_local - dot(a1_local, t) * t;
 a1_local = a1_local / norm(a1_local);

 a1(c,:) = a1_local; % store
 a2(c,:) = cross(t, a1_local);
end
end

function kappa = getkappa( q, m1, m2 )
% Inputs:
% q: DOF vector
% m1 is the first material director (size is ne x 3)
% m2 is the second material director (size is ne x 3)
% Output
% kappa is a matrix of size nv x 2. Each node has kappa_1 and kappa_2.
nv = ( length(q) + 1) / 4; % length(q) = 4*nv-1
ne = nv - 1;
kappa = zeros(nv, 2);
for c=2:ne % All internal nodes except first and last
 node0 = q(4*c-7:4*c-5);
 node1 = q(4*c-3:4*c-1);
 node2 = q(4*c+1:4*c+3);

 m1e = m1(c-1, :); % m1 vector of c-1-th edge
 m2e = m2(c-1, :); % m2 vector of c-1-th edge
 m1f = m1(c, :); % m1 vector of c-th edge
 m2f = m2(c, :); % m2 vector of c-th edge

 kappaLocal = computekappa(node0, node1, node2, m1e, m2e, m1f, m2f );

 kappa(c,1) = kappaLocal(1);
 kappa(c,2) = kappaLocal(2);
end
end

function d = parallel_transport(u, t1, t2)
% This function parallel transports a vector u from
% tangent t1 to t2.
b = cross(t1,t2);
if norm(b) == 0
 d = u;
else
 b = b / norm(b);
 % Good practice for safety
 b = b - dot(b,t1)*t1;
 b = b / norm(b);
 b = b - dot(b,t2)*t2;
 b = b / norm(b);

 n1 = cross(t1, b);
 n2 = cross(t2, b);
 d = dot(u,t1) * t2 + dot(u,n1) * n2 + ...
 dot(u,b) * b;
end
end

function [Fb, Jb] = getFb(q, m1, m2)
global kappaBar EI voronoiLength
nv = (length(q)+1) / 4;
Fb = zeros( size(q) );
Jb = zeros( length(q), length(q) );
for c=2:nv-1 % Compute bending force at each internel node

 node0 = transpose(q(4*c-7:4*c-5));
 node1 = transpose(q(4*c-3:4*c-1));
 node2 = transpose(q(4*c+1:4*c+3));

 m1e = m1(c-1, :); % m1 vector of c-1-th edge
 m2e = m2(c-1, :); % m2 vector of c-1-th edge
 m1f = m1(c, :); % m1 vector of c-th edge
 m2f = m2(c, :); % m2 vector of c-th edge

 [dF, dJ] = ...
 gradEb_hessEb(node0, node1, node2, ...
 m1e, m2e, m1f, m2f, ...
 kappaBar(c,:), voronoiLength(c), EI);

 ind = 4*c-7:4*c+3; % 11 numbers
 Fb(ind) = Fb(ind) - dF;
 Jb(ind, ind) = Jb(ind, ind) - dJ;

end
end

function [Ft, Jt] = getFt(q, refTwist)
global GJ voronoiLength
nv = (length(q)+1) / 4;
Ft = zeros( size(q) );
Jt = zeros( length(q), length(q) );
for c=2:nv-1 % Compute bending force at each internel node

 node0 = transpose(q(4*c-7:4*c-5));
 node1 = transpose(q(4*c-3:4*c-1));
 node2 = transpose(q(4*c+1:4*c+3));
 theta_e = q(4*c-4);
 theta_f = q(4*c);

 [dF, dJ] = ...
 gradEt_hessEt(node0, node1, node2, ...
 theta_e, theta_f, refTwist(c), ...
 voronoiLength(c), GJ);

 ind = 4*c-7:4*c+3; % 11 numbers
 Ft(ind) = Ft(ind) - dF;
 Jt(ind, ind) = Jt(ind, ind) - dJ;
end
end

function [Fs, Js] = getFs(q)
global EA refL
nv = (length(q)+1) / 4;
ne = nv - 1;
Fs = zeros( size(q) );
Js = zeros( length(q), length(q) );
for c=1:ne % Each edge

 node1 = transpose(q(4*c-3:4*c-1)); % c-th node
 node2 = transpose(q(4*c+1:4*c+3)); % c+1-th node

 [dF, dJ] = ...
 gradEs_hessEs(node1, node2, refL(c), EA);

 ind = [4*c-3, 4*c-2, 4*c-1, 4*c+1, 4*c+2,4*c+3]; % 6 numbers
 Fs(ind) = Fs(ind) - dF;
 Js(ind, ind) = Js(ind, ind) - dJ;

end
end

function plotrod(q, a1, a2, m1, m2, ctime)
nv = (length(q)+1)/4;
x1 = q(1:4:end);
x2 = q(2:4:end);
x3 = q(3:4:end);
L = sum(sqrt( (x1(2:end) - x1(1:end-1)).^2 + ...
 (x2(2:end) - x2(1:end-1)).^2 + ...
 (x3(2:end) - x3(1:end-1)).^2)) / 3;
a1 = 0.1*L * a1;
a2 = 0.1*L * a2;
m1 = 0.1*L * m1;
m2 = 0.1*L * m2;
h1 = figure(1);
% set(h1, 'visible', 'off');
clf()
plot3(x1,x2,x3, 'ko-');
hold on
plot3(x1(1),x2(1),x3(1), 'r^');
for c=1:nv-1
 xa = q(4*c-3:4*c-1);
 xb = q(4*c+1:4*c+3);
 xp = (xa+xb)/2;
 p1 = plot3( [xp(1), xp(1) + a1(c,1)], [xp(2), xp(2) + a1(c,2)], ...
 [xp(3), xp(3) + a1(c,3)], 'b--', 'LineWidth', 2);
 p2 = plot3( [xp(1), xp(1) + a2(c,1)], [xp(2), xp(2) + a2(c,2)], ...
 [xp(3), xp(3) + a2(c,3)], 'c--', 'LineWidth', 2);
 p3 = plot3( [xp(1), xp(1) + m1(c,1)], [xp(2), xp(2) + m1(c,2)], ...
 [xp(3), xp(3) + m1(c,3)], 'r-');
 p4 = plot3( [xp(1), xp(1) + m2(c,1)], [xp(2), xp(2) + m2(c,2)], ...
 [xp(3), xp(3) + m2(c,3)], 'g-');
end
hold off
legend([p1,p2,p3,p4], 'a_1','a_2','m_1','m_2');
title(num2str(ctime, 't=%f'));
axis equal
view(3);
xlabel('x');
ylabel('y');
zlabel('z');
drawnow
end