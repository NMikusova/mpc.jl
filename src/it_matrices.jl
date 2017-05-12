function it_matrices(A,B,C,D,x0,um1,Qy,Qu,ref,N,umin,umax,dumin,dumax,ymin,ymax)

## MODEL SIZE CHECK
modelt_size_check(A,B,C,D,x0,um1)

## NUMBER OF STATES
nx = size(A,1)
## NUMBER OF INPUTS
nu = size(B,2)
## NUMBER OF OUTPUTS
ny = size(C,1)

## Y = Ct*x0 + Dt*U
Ct = C
for i = 1:(N-1)
  Ct = [Ct; C*(A^i)]
end

c = D
for i = 0:(N-2)
  c = [c; C*(A^i)*B]
end
Dt = c
for i = 1:(N-1)
  ci = [zeros(ny*i,nu);c[1:(end-ny*i),:]]
  Dt = [Dt ci]
end

## dU = L*U + lambda*u-1
lambda = [-eye(nu,nu);zeros(nu*(N-1),nu)]

c = [eye(nu,nu);-eye(nu,nu);zeros(nu*(N-2),nu)]
L = c
for i = 1:(N-1)
  if i == (N-1)
    L = [L [zeros(nu*(N-1),nu);eye(nu,nu)]]
  else
    L = [L circshift(c, i*nu)]
  end
end

## min (Y-R)'*Qyt*(Y-R) + dU'*Qut*dU

weight_size_check(Qy,Qu,C,B)

Qyt = kron(eye(N),Qy)
Qut = kron(eye(N),Qu)

R = []
if length(ref) == ny
  for i = 1:N
    R = [R; ref]
  end
else
  error("Number of rows in C must be equal to number of rows in ref.")
end

## P,q,r
P = 2*(Dt'*Qyt*Dt + L'*Qut*L)
q = 2*(x0'*Ct'*Qyt*Dt - R'*Qyt*Dt + um1'*lambda'*Qut*L)'
r = x0'*Ct'*Qyt*Ct*x0 + R'*Qyt*R + um1'*lambda'*Qut*lambda*um1 - 2*x0'*Ct'*Qyt*R

## constraints

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

bieq = [Umax;-Umin]
Aieq = [eye(nu*N);-eye(nu*N)]

if isempty(ymin) == false && isempty(ymax) == false
  Ymax = [];Ymin = [];
  if length(ymin) == ny && length(ymax) == ny
    for i = 1:N
      Ymax = [Ymax; ymax]
      Ymin = [Ymin; ymin]
    end
  else
    error("Number of rows in C must be equal to number of rows in ymin and ymax.")
  end
  Aieq = [Aieq; Dt;-Dt]
  bieq = [bieq; Ymax-Ct*x0;-Ymin+Ct*x0]
end

if isempty(dumin) == false && isempty(dumax) == false
  dUmax = [];dUmin = [];
  if length(dumin) == nu && length(dumax) == nu
    for i = 1:N
      dUmax = [dUmax; dumax]
      dUmin = [dUmin; dumin]
    end
  else
    error("Number of columns in B must be equal to number of rows in dumin and dumax.")
  end
  Aieq = [Aieq; L;-L]
  bieq = [bieq; dUmax - lambda*um1;-dUmin + lambda*um1]
end

return P,q,r,Aieq,bieq

end
