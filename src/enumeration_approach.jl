function enumeration_approach(H,q,fu,theta,G,E,w)

output = []
m = size(G,1)
inactive1 = []
inactive = Array{Int}[]
itr = 0

## level 0
print("Level 0: Creating region...\n")
aU = -inv(H)*q'
#bU = zeros(size(aU,1),1)
bU = -inv(H)*fu
Ar = G*aU - E
br = w - G*bU
push!(output, region(Ar,br,aU,bU))

## level 1
print("Level 1: Processing...\n")
for i = 1:m
  itr = itr+1
  GA = G[[i],:]
  rN = collect(1:m)
  EA = E[[i],:]; wA = w[[i],:]
  deleteat!(rN, i)
  GN = G[rN,:]; EN = E[rN,:]; wN = w[rN,:]
  c,Aeq,beq,A,b = LPform(H,q,fu,theta,GA,EA,wA,GN,EN,wN,"s")
  cN,AN,bN = LPnullspace(c,Aeq,beq,A,b)
  sol = linprog(cN[1,:],AN,'<',bN[:,1],-Inf,Inf,ClpSolver())
  if sol.status == :Optimal
    print("Combination $itr is active: Creating region...\n")
    Ar,br,aU,bU = regions(H,q,fu,GA,EA,wA,GN,EN,wN)
    push!(output, region(Ar,br,aU,bU))
  else
    print("Combination $itr is inactive: ")
    c,Aeq,beq,A,b = LPform(H,q,fu,theta,GA,EA,wA,GN,EN,wN,"ws")
    cN,AN,bN = LPnullspace(c,Aeq,beq,A,b)
    sol = linprog(cN[1,:],AN,'<',bN[:,1],-Inf,Inf,ClpSolver())
    if sol.status == :Infeasible
      print("Removing all combinations containing current one...\n")
      inactive1 = [inactive1; i]
    else
      print("Removing combination...\n")
    end
  end
end

v = collect(1:m)
deleteat!(v, inactive1)

## levels 2 - m
for i = 2:m
  print("Level $i: Processing...\n")
  for j in combinations(v,i)
    if i>2
      test = false
      for k = 1:length(inactive)
        f = findin(j,inactive[k])
        if length(f) == length(inactive[k])
          test = true
          break
        end
      end
      if test == true
        continue
      end
    end
    itr = itr+1
    rA = j
    rN = collect(1:m)
    GA = G[rA,:]
    if rank(GA) == size(GA,1)
      EA = E[rA,:]; wA = w[rA,:]
      deleteat!(rN, rA)
      GN = G[rN,:]; EN = E[rN,:]; wN = w[rN,:]
      c,Aeq,beq,A,b = LPform(H,q,fu,theta,GA,EA,wA,GN,EN,wN,"s")
      cN,AN,bN = LPnullspace(c,Aeq,beq,A,b)
      sol = linprog(cN[1,:],AN,'<',bN[:,1],-Inf,Inf,ClpSolver())
      if sol.status == :Optimal
        print("Combination $itr is active: Creating region...\n")
        Ar,br,aU,bU = regions(H,q,fu,GA,EA,wA,GN,EN,wN)
        push!(output, region(Ar,br,aU,bU))
      else
        print("Combination $itr is inactive: ")
        c,Aeq,beq,A,b = LPform(H,q,fu,theta,GA,EA,wA,GN,EN,wN,"ws")
        cN,AN,bN = LPnullspace(c,Aeq,beq,A,b)
        sol = linprog(cN[1,:],AN,'<',bN[:,1],-Inf,Inf,ClpSolver())
        if sol.status == :Infeasible
          print("Removing all combinations containing current one...\n")
          push!(inactive, rA)
        else
          print("Removing combination...\n")
        end
      end
    else
      print("Combination $itr doesn't have full rank: Continue...\n")
    end
  end
end

print("Done!")
return output

end
