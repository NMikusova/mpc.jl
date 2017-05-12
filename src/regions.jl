function regions(H,q,fu,GA,EA,wA,GN,EN,wN)
  #regions(H,q,GA,EA,wA,GN,EN,wN,Y,theta)

alfaL = inv(GA*inv(H)*GA')*(-GA*inv(H)*q' - EA)
betaL = -inv(GA*inv(H)*GA')*(wA+GA*inv(H)*fu)

alfaU = inv(H)*(-q' - GA'*alfaL)
betaU = -inv(H)*fu -inv(H)*GA'*betaL

A = [GN*alfaU - EN; -alfaL]
b = [wN - GN*betaU; betaL]

# alfaJ = 1/2*alfaU'*H*alfaU + q*alfaU + Y
# betaJ = betaU'*H*alfaU + betaU'*q'
# gammaJ = 1/2*betaU'*H*betaU
#
# obj = theta'*alfaJ*theta + betaJ*theta + gammaJ

return A,b,alfaU,betaU

end
