## Regulation Example
## 2 states, 1 input

A = [-0.0490 -0.3278;0.2622 0.7376];
B = [0.2622; 0.2099];
Qx = eye(2,2);
Qu = 0.1;
N = 5;
x0 = [-1;3];
umin = -1.5;
umax = 0.5;
t = 20

tic()
mdl = rmodel(A,B,x0,Qx,Qu,umin,umax,[],[],N)
u,x,obj = iregulation(mdl,t)
toc()
