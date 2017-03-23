function rho = cossim(ws,w0)
% Calculate the cosine similarity between two vectors
rho = (ws'*w0)/(norm(ws)*norm(w0));
end