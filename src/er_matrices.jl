function er_matrices(A,B,x0,Qx,Qu,N,umin,umax,xmin,xmax)

  ## MODEL SIZE CHECK
  modelr_size_check(A,B,x0)
  weight_size_check(Qx,Qu,A,B)

  ## NUMBER OF STATES
  nx = size(A,1)
  ## NUMBER OF INPUTS
  nu = size(B,2)

  Qxt = kron(eye(N),Qx)
  Qut = kron(eye(N),Qu)

  At = eye(nx,nx)
  for i = 1:(N-1)
    At = [At; A^i]
  end

  c = zeros(nx,nu)
  for i = 0:(N-2)
    c = [c; (A^i)*B]
  end
  Bt = c
  for i = 1:(N-1)
    ci = [zeros(nx*i,nu);c[1:(end-nx*i),:]]
    Bt = [Bt ci]
  end

  H = 2*(Bt'*Qxt*Bt + Qut)
  theta = x0
  q = 2*At'*Qxt*Bt
  Y = At'*Qxt*At

  if isempty(umin) == true || isempty(umax) == true
    error("umin and umax must be defined.")
  end

  Umin = Array{Float64}(0,size(B,2))
  Umax = Array{Float64}(0,size(B,2))
  if length(umin) == nu && length(umax) == nu
    for i = 1:N
      Umax = [Umax; umax]
      Umin = [Umin; umin]
    end
  else
    error("Number of columns in B must be equal to number of rows in umin and umax.")
  end

  if isempty(xmin) == false && isempty(xmax) == false
    Xmax = []; Xmin = [];
    if length(xmin) == nx && length(xmax) == nx
      for i = 1:N
        Xmax = [Xmax; xmax]
        Xmin = [Xmin; xmin]
      end
    else
      error("Number of rows in A must be equal to number of rows in xmin and xmax.")
    end
    G = [eye(N*nu);-eye(N*nu);Bt;-Bt]
    E = [zeros(2*N*nu,nx);-At;At]
    w = [Umax;-Umin;Xmax;-Xmin]
  else
    G = [eye(N*nu);-eye(N*nu)]
    E = zeros(2*N*nu,nx)
    w = [Umax;-Umin]
  end

  return H,q,Y,theta,G,E,w

end
