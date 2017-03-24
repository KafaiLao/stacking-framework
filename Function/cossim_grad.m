function rho = cossim_grad(ws,w0)
rho = w0/(norm(ws)*norm(w0))-(ws*ws'*w0)/(norm(ws)^3*norm(w0));