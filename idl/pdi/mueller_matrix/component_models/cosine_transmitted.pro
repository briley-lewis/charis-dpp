FUNCTION cosine_transmitted, th_i, n1, n2
cos_t = sqrt(1.-((n1/n2)*sin(th_i))^2)
return, cos_t
end