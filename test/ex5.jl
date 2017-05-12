## Tracking Example
## 2 states, 1 output, 1 input

A = [-0.0490 -0.3278;0.2622 0.7376];
B = [0.2895;0.2284];
C = [1 1.5];
D = [0];
x0 = [0;0];
um1 = 0;
N = 5;
ref = 2;
Qy = 1e3;
Qu = 1e-2;
umin = -2.2
umax = 2.2
dumin = -0.75
dumax = 1
t = 20

tic()
mdl = tmodel(A,B,C,D,x0,um1,Qy,Qu,ref,umin,umax,dumin,dumax,[],[],N)
u,y,x,r = itracking(mdl,t)
toc()
