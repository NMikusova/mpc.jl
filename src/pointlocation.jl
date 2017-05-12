function pointlocation(r,theta,nu)

u = collect(1:nu)*0

for i = 1:length(r)
  if maximum(r[i].A*theta - r[i].b) <= 1e-12
    u = (r[i].aU*theta + r[i].bU)[1:nu]
    break
  end
end

return u

end
