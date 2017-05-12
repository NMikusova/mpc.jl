function ir_matrices(A,B,x0,Qx,Qu,N,umin,umax,xmin,xmax)

  ## MODEL SIZE CHECK
  modelr_size_check(A,B,x0)

  ## NUMBER OF STATES
  nx = size(A,1)
  ## NUMBER OF INPUTS
  nu = size(B,2)

  weight_size_check(Qx,Qu,A,B)

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
  q = (2*x0'*At'*Qxt*Bt)'
  r = x0'*At'*Qxt*At*x0

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
    A = [eye(N*nu);-eye(N*nu);Bt;-Bt]
    b = [Umax;-Umin;Xmax-At*x0;-Xmin+At*x0]
  else
    A = [eye(N*nu);-eye(N*nu)]
    b = [Umax;-Umin]
  end

  return H,q,r,A,b

end
