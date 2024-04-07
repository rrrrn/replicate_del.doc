fpprime_theta <- function(p, t, a, u, Omega, vec, g_EL){
  temp = 1-t[1]*(x>1.517427) + t[1]*(rep(p, length(x)))-t[2]*y+t[2]*rep(2, length(y))-t[3]*z+t[3]*x
  temp = sum(1/(temp**2))
  return(temp + 2*Omega)
}