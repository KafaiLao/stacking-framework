function [Wall_new,v_new,w0_new,ls_new,rho_new,dws,dv,dw0] = updateGrad_fixV(Zall,Y,Wall,v,w0,beta,alpha,mu,ls,rho)
%Zall [Sample x channel x subejct]
%Y [Sample x number of harmonic*2]
%Wall [spatial filter size x subejct]
%Vall [number of harmonic*2 x subject]
%w0 [spatial filter size x 1]
N = size(Zall,1)-1;
S = size(Zall,3);
dw0 = zeros(size(w0));
dws = zeros(size(Wall));
dv = zeros(size(v));
dls = zeros(size(ls));

for subject = 1:S
    Zs = squeeze(Zall(:,:,subject));
    ws = Wall(:,subject);
    
    comp = w0/(norm(ws)*norm(w0)) - (ws*ws'*w0)/(norm(ws)^3*norm(w0));
    dws(:,subject) = (1-beta)/(N*S)*Zs'*Y*v + beta/S*comp - 2*mu*(ws'*(Zs'*Zs)*ws/N - 1)*(Zs'*Zs)*ws/N + 2*ls(subject)*(Zs'*Zs)*ws/N;
    dv = dv + (1-beta)/(N*S)*Y'*Zs*ws; 
    comps = ws/(norm(ws)*norm(w0)) - (w0*w0'*ws)/(norm(w0)^3*norm(ws));
    dw0 = dw0 + beta/S*comps;
    dls(subject) = -mu*(ws'*(Zs'*Zs)*ws/N - 1);
end
dv = dv - 2*mu*(v'*(Y'*Y)*v/N - 1)*(Y'*Y)*v/N + 2*rho*(Y'*Y)*v/N;
drho = -mu*(v'*(Y'*Y)*v/N - 1);

Wall_new = Wall + alpha*dws;
v_new = v + alpha*dv;
w0_new = w0 + alpha*dw0;
ls_new = ls + dls;
rho_new = rho + drho;


