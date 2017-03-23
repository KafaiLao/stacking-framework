function output = objFun_fixV(Zall,Y,Wall,v,w0,beta)
%Zall [Sample x channel x subejct]
%Y [Sample x number of harmonic*2]
%Wall [spatial filter size x subejct]
%V [number of harmonic*2 x 1]
%w0 [spatial filter size x 1]
output = 0;
S = size(Zall,3);
for subject = 1:S
    Zs = squeeze(Zall(:,:,subject));
    ws = Wall(:,subject);
    corr_num = ws'*Zs'*Y*v;
    corr_den = sqrt(ws'*(Zs'*Zs)*ws)*sqrt(v'*(Y'*Y)*v);
    output = output + (1-beta)/S*corr_num/corr_den + beta/S*cossim(ws,w0);
end