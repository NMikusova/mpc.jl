function LPform(H,q,fu,x0,GA,EA,wA,GN,EN,wN,stat)

## TOTAL number of inputs
nu = size(H,1)
## number of states
nx = length(x0)
## number of active constraints
cA = size(GA,1)

## original problem:
## min -t
## s.t. H*U + q'*tht + GA'*lA + fu = 0 !optional
##      GA*U - EA*tht = wA
##      t + GN*U - EN*tht <= wN
##      t - lA <= 0
##      -t <= 0

## standard LP form:
## min c'*x
## s.t. Aeq*x = beq
##      A*x <= b

## x = [t;U;tht;lA]

## objective function
c = [-1; zeros(nu,1); zeros(nx,1); zeros(cA,1)]

# inequalities
cN = size(GN,1) # number of inactive constraints
A = [zeros(cN,1) GN -EN zeros(cN,cA); ones(cA) zeros(cA,nu) zeros(cA,nx) -eye(cA); -1 zeros(1,nu) zeros(1,nx) zeros(1,cA)]
b = [wN; zeros(cA,1); 0]

# equalities
if stat == "s"
  Aeq = [zeros(nu,1) H q' GA'; zeros(cA,1) GA -EA zeros(cA,cA)]
  beq = [fu;wA]
else
  Aeq = [zeros(cA,1) GA -EA zeros(cA,cA)]
  beq = wA
end

return c,Aeq,beq,A,b

end
