l = .155;
g = 9.81;
rho =11;
Ts=0.01;
A = [1,0,0;
    0,1,Ts/l;
    0,-g*Ts,1];
B = [Ts;
    -Ts/l;
    0];
Cz = [1, l, 0]; % For output is position tip of pendulum LQR
syms z
Q = rho*transpose(Cz)*Cz;
R=1
K_FB_31=dlqr(A,B,Q,R)

C =[1 0 0;0 -1 0];

H = C*inv(z*eye(3)-A+B*K_FB_31)*B;
H_y1 = H(1);
H_y2 = H(2);

lim=limit(H_y1,z,1);
double(lim)
K_pre = 1/double(lim)