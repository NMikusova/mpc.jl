# mpc.jl
Julia package for Implicit and Explicit Model Predictive Control

## Download and Install
In console: `Pkg.clone("git://github.com/NMikusova/mpc.jl.git")`

mpc.jl depends on these registered packages (installed automatically):
* MathProgBase
* Clp
* Ipopt
* Combinatorics


## How to use the tool
Initialization:
`using mpc`

### Regulation Problem

Objective: regulate the system toward its origin, i.e., to a zero state.

Input arguments:

Argument | Description
------------ | -------------
`A,B` | State-space matrices
`x0` | Initial condition for states
`Qx,Qu` | Weighting matrices
`umin,umax` | Box constraints on inputs
`xmin,xmax` | Box constraints on states
`N` | Prediction horizon
`t` | Simulation time

Output arguments:

Argument | Description
------------ | -------------
`u` | Optimal control action
`x` | Optimal states

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

### Tracking Problem

Objective: offset-free tracking of non-zero reference

Input arguments:

Argument | Description
------------ | -------------
`A,B,C,D` | State-space matrices
`x0` | Initial condition for states
`um1` | Initial condition for inputs
`Qy,Qu` | Weighting matrices
`ref` | Reference
`umin,umax` | Box constraints on inputs
`dumin,dumax` | Box constraints on input increments
`ymin,ymax` | Box constraints on outputs
`N` | Prediction horizon
`t` | Simulation time

Output arguments:

Argument | Description
------------ | -------------
`u` | Optimal control action
`y` | Optimal outputs
`x` | Optimal states

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
