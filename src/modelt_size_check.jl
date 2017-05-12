function modelt_size_check(A,B,C,D,x0,um1)

modelr_size_check(A,B,x0)

if size(A,2) != size(C,2)
  error("A and C must be matrices with the same number of columns.")
end

if size(B,2) != size(D,2)
  error("B and D must be matrices with the same number of columns.")
end

if size(C,1) != size(D,1)
  error("C and D must be matrices with the same number of rows.")
end

if size(B,2) != size(um1,1)
  error("Number of columns in B must be equal to number of rows in um1.")
end

end
