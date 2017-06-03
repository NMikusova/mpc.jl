# mpc.jl
Julia package for Implicit and Explicit Model Predictive Control

## Download and Install
In console: `Pkg.clone("git://github.com/NMikusova/mpc.jl.git")`

## How to use the tool
Initialization:
`using mpc`

### Regulation Problem

Input arguments:

Argument | Description
------------ | -------------
A,B | State-space matrices
x0 | Initial condition
Qx,Qu | Weighting matrices
umin,umax | Box constraints on inputs
xmin,xmax | Box constraints on states
N | Prediction horizon

#### Regulation problem + Implicit MPC
```
mdl = rmodel(A,B,x0,Qx,Qu,umin,umax,xmin,xmax,N)
u,x = iregulation(mdl,t)
```

`xmin` and `xmax` are optional:
```
mdl = rmodel(A,B,x0,Qx,Qu,umin,umax,[],[],N)
```

#### Regulation problem + Explicit MPC
```
mdl = rmodel(A,B,x0,Qx,Qu,umin,umax,xmin,xmax,N)
u,x = eregulation(mdl,t)
```

`xmin` and `xmax` are optional:
```
mdl = rmodel(A,B,x0,Qx,Qu,umin,umax,[],[],N)
```

#### Tracking problem + Implicit MPC
```
mdl = tmodel(A,B,C,D,x0,um1,Qy,Qu,ref,umin,umax,dumin,dumax,ymin,ymax,N)
u,y,x = itracking(mdl,t)
```

`dumin`,`dumax`,`ymin` and `ymax` are optional:
```
mdl = tmodel(A,B,C,D,x0,um1,Qy,Qu,ref,umin,umax,[],[],[],[],N)
```

#### Tracking problem + Explicit MPC
```
mdl = tmodel(A,B,C,D,x0,um1,Qy,Qu,ref,umin,umax,dumin,dumax,ymin,ymax,N)
u,y,x = etracking(mdl,t)
```

`dumin`,`dumax`,`ymin` and `ymax` are optional:
```
mdl = tmodel(A,B,C,D,x0,um1,Qy,Qu,ref,umin,umax,[],[],[],[],N)
```
