function weight_size_check(Qxy,Qu,AC,B)

## WEIGHT MATRICES SIZE CHECK
if size(Qxy,1) != size(Qxy,2)
  error("Qx/Qy must be square matrix.")
end

if size(Qu,1) != size(Qu,2)
  error("Qu must be square matrix.")
end

if size(AC,1) != size(Qxy,1)
  error("A and Qx/C and Qy must be matrices with the same number of rows.")
end

if size(B,2) != size(Qu,1)
  error("Number of columns in B must be equal to number of rows in Qu.")
end

end
