function output = objFun_matlab(coef,Zall,Y,beta)
%Zall [Sample x channel x subejct]
%Y [Sample x number of harmonic*2]
%Wall [spatial filter size x subejct]
%Vall [number of harmonic*2 x subject]
%w0 [spatial filter size x 1]
temp = coef(1:9*35);
Wall = reshape(temp,9,35);
temp = coef(9*35+1:9*35+10*35);
Vall = reshape(temp,10,35);
w0 = coef(9*35+10*35+1:end);

output = 0;
S = size(Zall,3);
for subject = 1:S
    Zs = squeeze(Zall(:,:,subject));
    ws = Wall(:,subject);
    vs = squeeze(Vall(:,subject));
    fnum = ws'*Zs'*Y*vs;
    fden = sqrt(ws'*(Zs'*Zs)*ws)*sqrt(vs'*(Y'*Y)*vs);
    output = output + fnum/fden + beta*cossim(ws,w0);
end
output = -output;