function iregulation(mdl,t)

## initialization
u = zeros(size(mdl.B,2),t)
x = zeros(size(mdl.A,1),t+1)
x[:,1] = mdl.x0
obj = zeros(t,1)

for k = 1:t
  P,q,r,Ai,bi = ir_matrices(mdl.A,mdl.B,x[:,k],mdl.Qx,mdl.Qu,mdl.N,mdl.umin,mdl.umax,mdl.xmin,mdl.xmax)
  ## activeset
  sol = quadprog(q[:,1],P,Ai,'<',bi[:,1],-Inf,Inf,IpoptSolver())
  obj[k] = sol.objval
  u[:,k] = sol.sol[1:size(mdl.B,2)]
  x[:,k+1] = mdl.A*x[:,k] + mdl.B''*u[:,k]
end

return u,x,obj

end
