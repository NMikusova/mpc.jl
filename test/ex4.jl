## Regulation Example
## 1 state, 1 input

A = 0.8;
B = 1;
x0 = 5;
N = 4;
Qx = 1;
Qu = 0.1;
umin = -1;
umax = 1;
t = 20

tic()
mdl = rmodel(A,B,x0,Qx,Qu,umin,umax,[],[],N)
u,x,r = eregulation(mdl,t)
toc()
