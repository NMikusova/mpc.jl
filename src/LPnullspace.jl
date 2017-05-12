function LPnullspace(c,Aeq,beq,A,b)

## standard LP form:
## min c'*x
## s.t. Aeq*x = beq
##      A*x <= b

## Aeq*x0 = beq
## A*F = 0
## x0 = pinv(Aeq)*beq
## F = nullspace(Aeq)
## x = x0 + F*alfa

## LP without equalities
## min c'*F*alfa + c'*x0
## s.t. A*F*alfa <= b - A*x0

x0 = pinv(Aeq)*beq
F = nullspace(Aeq)
cNew = c'*F
ANew = A*F
bNew = b - A*x0

return cNew,ANew,bNew

end
