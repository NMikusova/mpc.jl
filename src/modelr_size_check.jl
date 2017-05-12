function modelr_size_check(A,B,x0)

## MODEL SIZE CHECK
if size(A,1) != size(A,2)
  error("A must be square matrix.")
end

if size(A,1) != size(B,1)
  error("A and B must be matrices with the same number of rows.")
end

if size(A,1) != size(x0,1)
  error("A and x0 must be matrices with the same number of rows.")
end

end
