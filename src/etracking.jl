function etracking(mdl,t)

  H,q,theta,fu,Y,f0,f,G,E,w = et_matrices(mdl.A,mdl.B,mdl.C,mdl.D,mdl.x0,mdl.um1,mdl.Qy,mdl.Qu,mdl.ref,mdl.N,mdl.umin,mdl.umax,mdl.dumin,mdl.dumax,mdl.ymin,mdl.ymax)
  r = enumeration_approach(H,q,fu,theta,G,E,w)

  u = zeros(size(mdl.B,2),t+1)
  y = zeros(size(mdl.C,1),t)
  x = zeros(size(mdl.A,1),t+1)
  nx = size(mdl.x0,1)
  nu = size(mdl.um1,1)
  u[:,1] = mdl.um1
  x[:,1] = theta[1:nx]
  #
  for k = 1:t
    u[:,k+1] = pointlocation(r,theta,nu)
    y[:,k] = mdl.C*x[:,k] + mdl.D''*u[:,k+1]
    x[:,k+1] = mdl.A*x[:,k] + mdl.B''*u[:,k+1]
    theta = [x[:,k+1];u[:,k+1]]
  end

  return u,y,x,r

end
