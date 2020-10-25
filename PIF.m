function [res,c_est,b_est_mean,sn_est,b_est_mode,W,lkl] = PIF(z,aord,init,lam,alpha,beta,niter,Nbur)
%% low level parameters
sc=1;
N=length(z);
a=zeros(N,1);
b=binornd(2,lam,[N,1]);
b=ones(size(b));
s =[init; z];s11=s;
S=zeros(N,aord);
for j=aord:length(s)-1
    S(j-aord+1,:)=s(j:-1:j-aord+1)';
end
s = s(aord+1:end);
c=-pinv(S)*s;
%%
sn=1e-2;

best=zeros(N,niter);
aest=zeros(N,niter);
snest=zeros(1,niter);
cest=zeros(length(c),niter);

ep=1e-3;
wsz=3;nz=0;

win=[ep*ones(nz,1);rectwin(wsz)-ep];
A1=convmtx(win,length(z));%A=A1;
A = A1(1:end-wsz+1,:);%A=A';
ee = ones(1,length(z));
d =s+S*c;
ep=1e-5;
for k=1:niter
    D=diag(d.^2);
    z = A1*D*ee'/(2*sn);
    z=z(wsz:end);
    j=1;
    while(j<=length(z))
        b(j)=0;
        tmp33=(sum(log(A*b+ep)));
        b(j)=1;
        t1=sqrt(exp(-sum(log(A*b+ep))+tmp33));
        
        p1(j)=exp(-z(j));
        ttt=(1-lam)*t1;
        p(j)=p1(j)*lam/(p1(j)*lam+ttt);
        if(rand(1)>p(j))
            b(j)=0;
        end
        j=j+1;
    end
    w=A*b+ep;
    W=diag(w);
    W1=W/sn;
    SS=S'*W1*S;
    SSS=(SS+(eye(size(S,2))/sc));
    Sc = (inv(SSS));Sc=(Sc+Sc')/2;
    mu = -Sc*S'*W1*s;
    c = mvnrnd(mu,Sc)';
    d =s+S*c;
    val =d'*W*d;

    sn = 1./gamrnd((N/2)+alpha,1./((.5*(val))+beta));

    best(:,k)=b;
    aest(:,k)=a;
    cest(:,k)=c;
    west(:,k)=w;
    snest(k)=(sn);
    lkl(k)=-(S*c+s)'*(W1/2)*(S*c+s)-sum(log(diag(W)+1e-40))/2;%-(1/2)*log((2*pi)^length(sn) * sn);
end
c_est=[1;mean(cest(:,Nbur:end),2)];
b_est_mean = detectGCI(best,Nbur,'mean');
b_est_mode = detectGCI(best,Nbur,'mode');
sn_est=mean(snest(Nbur:end));
res=s11-filter([0;-c_est(2:end)],1,s11);
res=res(aord+1:end);
W=w;
% res=filter(1,[1,-0.98],res);
% res=res./min(res);
end


function b_est=detectGCI(best,Nbur,mtd)

switch mtd
    case 'mean'
        b_est = mean(best(:,Nbur:end),2);
    case 'mode'
        b_est=mode(best,2);
end
end