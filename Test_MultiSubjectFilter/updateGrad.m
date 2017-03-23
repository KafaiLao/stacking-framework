function [Wall_new,Vall_new,w0_new,ls_new,ps_new,dws,dvs,dw0] = updateGrad(Zall,Y,Wall,Vall,w0,beta,alpha,mu,ls,ps)
%Zall [Sample x channel x subejct]
%Y [Sample x number of harmonic*2]
%Wall [spatial filter size x subejct]
%Vall [number of harmonic*2 x subject]
%w0 [spatial filter size x 1]
N = size(Zall,1)-1;
S = size(Zall,3);
dw0 = zeros(size(w0));
dws = zeros(size(Wall));
dvs = zeros(size(Vall));
dls = zeros(size(ls));
dps = zeros(size(ps));
for subject = 1:S
    Zs = squeeze(Zall(:,:,subject));
    ws = Wall(:,subject);
    vs = squeeze(Vall(:,subject));
    
    comp = w0/(norm(ws)*norm(w0)) - (ws*ws'*w0)/(norm(ws)^3*norm(w0));
    dws(:,subject) = (1-beta)/(N*S)*Zs'*Y*vs + beta/S*comp - 2*mu*(ws'*(Zs'*Zs)*ws/N - 1)*(Zs'*Zs)*ws/N + 2*ls(subject)*(Zs'*Zs)*ws/N;
    dvs(:,subject) = (1-beta)/(N*S)*Y'*Zs*ws - 2*mu*(vs'*(Y'*Y)*vs/N - 1)*(Y'*Y)*vs/N + 2*ps(subject)*(Y'*Y)*vs/N;
    comps = ws/(norm(ws)*norm(w0)) - (w0*w0'*ws)/(norm(w0)^3*norm(ws));
    dw0 = dw0 + beta/S*comps;
    dls(subject) = -mu*(ws'*(Zs'*Zs)*ws/N - 1);
    dps(subject) = -mu*(vs'*(Y'*Y)*vs/N - 1);
end

Wall_new = Wall + alpha*dws;
Vall_new = Vall + alpha*dvs;
w0_new = w0 + alpha*dw0;
ls_new = ls + dls;
ps_new = ps + dps;


