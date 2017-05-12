function itracking(mdl,t)

## initialization
u = zeros(size(mdl.B,2),t+1)
y = zeros(size(mdl.C,1),t)
x = zeros(size(mdl.A,1),t+1)
x[:,1] = mdl.x0
u[:,1] = mdl.um1
obj = zeros(t,1)

for k = 1:t
  P,q,r,Ai,bi = it_matrices(mdl.A,mdl.B,mdl.C,mdl.D,x[:,k],u[:,k],mdl.Qy,mdl.Qu,mdl.ref,mdl.N,mdl.umin,mdl.umax,mdl.dumin,mdl.dumax,mdl.ymin,mdl.ymax)
  ## activeset
  sol = quadprog(q[:,1],P,Ai,'<',bi[:,1],-Inf,Inf,IpoptSolver())
  obj[k] = sol.objval
  u[:,k+1] = sol.sol[1:size(mdl.B,2)]
  y[:,k] = mdl.C*x[:,k] + mdl.D''*u[:,k+1]
  x[:,k+1] = mdl.A*x[:,k] + mdl.B''*u[:,k+1]
end

return u,y,x,obj

end
