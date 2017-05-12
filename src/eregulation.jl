function eregulation(mdl,t)

  H,q,Y,theta,G,E,w = er_matrices(mdl.A,mdl.B,mdl.x0,mdl.Qx,mdl.Qu,mdl.N,mdl.umin,mdl.umax,mdl.xmin,mdl.xmax)
  fu = zeros(size(H,1),1)
  r = enumeration_approach(H,q,fu,theta,G,E,w)

  u = zeros(size(mdl.B,2),t)
  x = zeros(size(mdl.A,1),t+1)
  nu = size(mdl.um1,1)
  x[:,1] = theta
  #
  for k = 1:t
    u[:,k] = pointlocation(r,theta,nu)
    theta = mdl.A*theta + mdl.B''*u[:,k]
    x[:,k+1] = theta
  end

  return u,x,r

end
