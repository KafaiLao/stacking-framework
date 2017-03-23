function filterBank = genFB(startFreq,endFreq,stopFreq,fsample,Nfilter)
% [Input]
% startFreq: minimum pass frequency (e.g., the minimum stimulus frequency)
% endFreq: end frequency (e.g., the maximum stimulus frequency)
% stopFreq: maximum frequency for the filter bank
% fsample: sampling frequency
% Nfilter: number of filter bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startFreq = floor(startFreq);
B = endFreq - startFreq;
filterBank = cell(Nfilter,1);
Rp = 0.5; % Passband ripple
Rs = 60; % Stopband attenuation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design chebyshev filter
for n = 1:Nfilter
    Wp = [startFreq + (n-1)*B stopFreq];
    Ws = [startFreq + (n-1)*B-1 stopFreq+1];
    if Wp(1)/Ws(1) > Ws(2)/Wp(2)
        Ws(1) = Wp(2)/Ws(2)*Wp(1);
    elseif Wp(1)/Ws(1) < Ws(2)/Wp(2)
        Ws(2) = Wp(1)/Ws(1)*Wp(2);
    end
    % Evaluate the order of prototype filter
    omega = (Ws(2) - Ws(1))/(Wp(2) - Wp(1));
    order = ceil(acosh(sqrt((10^(-0.1*Rs) - 1)/(10^(-0.1*Rp) - 1)))/acosh(omega));
    [z,p,k] = cheb1ap(order,abs(Rp));
    % Analog prototype (low-pass) filter
    [bp,ap] = zp2tf(z,p,k);
    % Prewarpping frequencies for bilinear transform
    Wp(1) = 2*fsample*tan(Wp(1)*(2*pi/fsample)/2);
    Wp(2) = 2*fsample*tan(Wp(2)*(2*pi/fsample)/2);
    w0 = sqrt(Wp(1)*Wp(2)); %Center frequency
    BW = Wp(2) - Wp(1); %Bandwidth
    % Design analog filter
    [bt,at] = lp2bp(bp,ap,w0,BW);
    % Bilinear transform
    [b,a] = bilinear(bt,at,fsample);
    filterBank{n} = [b;a];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Last modified 23/03/2017 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Ka Fai Lao, University of Macau %%%%%%%%%%%%%%%%%%%%%