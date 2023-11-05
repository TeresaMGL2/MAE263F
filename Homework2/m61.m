x1 = [0,0,0]';
x2 = [0.5,0,0]';
x3 = [0.75,0.25,0]';
x4 = [0.75,0.5,0.25]';
m1 = [0,0,1]';
m2 = [0,0,1]';
m3 = [0,-1/sqrt(2),1/sqrt(2)]';

t1 = tangent(x1,x2);
t2 = tangent(x2,x3);
t3 = tangent(x3,x4);

v1 = signed(para(m1,t1,t2),m2,t2);
v2 = signed(para(m2,t2,t3),m3,t3);

function v = signed(u,m,t)
w = cross(u,m); %calculate the amount of rotation(length corresponds to sin)
v = atan2(norm(w),dot(u,m)); %cross product w is the sin length while the dot between ref frame u and material frame m is the cos
if (dot(t,w)<0) %check direction of angle
    v = -v;
end
end

function up = para(u,t1,t2)
b = cross(t1,t2);
b = b/norm(b);
n1 = cross(t1,b);
n2 = cross(t2,b);
up = dot(u,t1)*t2+dot(u,n1)*n2+dot(u,b)*b;
end

function t = tangent(x1,x2)
t = x2-x1;%ek
t = t/norm(t);
end
