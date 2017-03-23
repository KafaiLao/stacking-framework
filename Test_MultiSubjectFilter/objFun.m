function output = objFun(Zall,Y,Wall,Vall,w0,beta)
%Zall [Sample x channel x subejct]
%Y [Sample x number of harmonic*2]
%Wall [spatial filter size x subejct]
%Vall [number of harmonic*2 x subject]
%w0 [spatial filter size x 1]
output = 0;
S = size(Zall,3);
for subject = 1:S
    Zs = squeeze(Zall(:,:,subject));
    ws = Wall(:,subject);
    vs = squeeze(Vall(:,subject));
    corr_num = ws'*Zs'*Y*vs;
    corr_den = sqrt(ws'*(Zs'*Zs)*ws)*sqrt(vs'*(Y'*Y)*vs);
    output = output + (1-beta)/S*corr_num/corr_den + beta/S*cossim(ws,w0);
end