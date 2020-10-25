clear;
rng('default');
[x,fs]=audioread('arctic_b0427.wav');
x=x(:,1);
x=resample(x,1,2);
fs=fs/2;
xorg=x;
%% hpf as per IAIF paper
% Fstop = 40;                 % Stopband Frequency
% Fpass = 70;                 % Passband Frequency
% Nfir = round(300/16000*fs); % FIR numerator order
% if mod(Nfir,2) == 1
%     Nfir = Nfir + 1;
% end
% B = hpfilter_fir(Fstop,Fpass,fs,Nfir);
%%
x=filter([1,-0.98],1,x); %%pre-emphasis
Nw=30e-3*fs;
Ns=10e-3*fs;
alpha=Nw;
beta=1e-4;

niter=100;
Nbur=20;
p=8;
% j=p+1;
j=21000; %just for test
gf=zeros(size(x));
f0_all=SRH_PitchTracking(xorg,fs,100,350); %%

while((j+Nw)<23000)
    %% PWLP
    z=x(j:j+Nw);
    f0=f0_all(round(j/Ns)+1);
    init=x(j-p:j-1);
     T0i=ceil(fs/f0);
    %setting the lambda function of f0 estimate (if not required set it to 1e-4)
    lam=3*floor(length(z)/T0i)/length(z);
    [~,c_est,b_est_mean,sn_est,b_est_mode,W,lkl]= PIF(z,p,init,lam,alpha,beta,niter,Nbur);
    gfd=filter(c_est,1,xorg(j-p:j+Nw));
    gfd=gfd(p+1:end);
    gf = filter(1,[1 -0.98],gfd(end:-1:1)); % Lip radiation compensation
    gf = -gf(end:-1:1);
    gfd=diff(gf);
    
    subplot(311);
    plot(xorg(j:j+Nw),'r');
    title('raw speech');
    subplot(312);
    plot(gfd,'r');hold on;plot(W*0.1,'m');hold off;
    title('flow derivative');
    subplot(313);
    plot(gf,'r');
    hold off;
    title('flow');
    legend('pwlp');
    
    j=j+Ns;
    
end
% plot(xorg);hold on;plot(filter(1,[1,-0.98],gf),'r');hold off;
