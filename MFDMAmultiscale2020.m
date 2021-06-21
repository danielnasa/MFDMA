clear all
clc;
close all

jj=0;
step=100; % 50 points are 1 days
for w=1:step:2901
    taux=w;
    www=w
    jj=jj+1
    
files=dir('PDCQ4withNoiseAndPlanets.txt')

s=length(files);

data=load(files.name,'-ascii');
    
    x=data(:,2); %DANGER see dimension of matrix
    %x=x-1; %to use when data aren't normalized in ZERO
    %x=x/nanmedian(x)-1;
    t=data(:,1);
    
    array=[t x];

array(any(isnan(array),2),:) = [];

x=array(:,2);
t=array(:,1);

%x=shuffle(x,200);
x=phaseran(x,200);

L=length(x);

    for j=1:L-taux
        x(j)=x(j+taux)-x(j);
    end

tauH(jj,1)=taux*29.4/(60*24);

nmin=10;
nmax=round(length(x)/10);
N=30;
nt=length(x);
q=-5:0.5:5;
theta=0;
%x=x(1:nt-1,1)-diff(x);
%x=diff(x);
%-------------

if size(x,2) == 1
x = x';
end
M = length(x);
MIN = log10(nmin) ;
MAX = log10(nmax) ;
n = (unique(round(logspace(MIN,MAX,N))))';
% Cons t ruct the cumulat i v e sum y
y = cumsum(x) ;
for i = 1:length(n)
lgth = n(i,1) ;
% Cal cu l at e the moving av e rag e func t ion \ wi d e t i l d e {y}
y1 = zeros(1,M-lgth+1) ;
for j = 1:M-lgth+1
y1(j) = mean(y(j:j+lgth-1)) ;
end
% Determine the r e s i d u a l e
eo=y(max(1,floor(lgth*(1-theta))):max(1,floor(lgth*(1-theta)))+length(y1)-1)-y1 ;
% Es t imate the root?mean?s quar e func t ion F10
for k=1:floor(length(eo)/lgth)
    F{i}(k)=sqrt(mean(eo((k-1)*lgth+1:k*lgth).^2));
end
end
% Cal cu l at e the q?th order o v e r a l l f l u c t u a t i o n func t ion Fq
for i=1:length(q)
for j =1:length(F)
fo=F{j};
if q(i)==0
Fq(j,i)=exp(0.5*mean(log(fo.^2 ))) ;
else
Fq(j,i)=(mean(fo.^q (i)))^(1/q(i));
end
end
end


% Cal cu l at e the mu l t i f r a c t a l s c a l i n g exponent tau ( q )
for i=1:size(Fq,2)
fq=Fq(:,i) ;
r=regstats(log(fq),log(n),'linear',{'tstat'}) ;
k=r.tstat.beta(2);
h(i,1)=k ;
end
tau=h.*q'-1;
hq=diff(tau)./(q(2)-q(1));
% Cal cu l at e the s i n g u l a r i t y s t r e n g t h func t ion alpha ( q ) and spectrum f ( alpha )
dx=7;
dx=fix((dx-1)/2);
for i=dx+1:length(tau)-dx
xx=q (i-dx:i+dx);
yy=tau(i-dx:i+dx);
r=regstats(yy,xx,'linear',{'tstat'}) ;
alpha(i,1)=r.tstat.beta(2);
end
alpha=alpha(dx+1:end) ;
fo=q(dx+1:end-dx)'.*alpha-tau(dx+1:end-dx);
Dq=(q(1:end-1)'.*hq)-tau(1:end-1);
fitalfa = polyfit(alpha,fo,3);
fittau = polyfit(q',tau,3);

%alphamaxO(w,1)=alpha(fo==1);
%maxalphaO=max(alpha(fo>=0));

%contriDiversityO(w,1)=maxalphaO-min(alpha); % degree of multifractality
%assymoriginalO(w,1)=(maxalphaO-alpha(fo==1))/(alpha(fo==1)-min(alpha)); % asymmetry
%dimDiversityLeftO(w,1)=max(fo)-min(fo(alpha==min(alpha))); % left side diversity
%dimDiversityRightO(w,1)=max(fo)-min(fo(alpha==maxalphaO)); % right side diversity
%CO(w,1)=min(fo(alpha==min(alpha)))-min(fo(alpha==maxalphaO));

for g=1:1:length(q)
    
    hurstO(g,jj)=h(g); % Hurst exponent
    
end
end
%MultifIndex=[CO alphamaxO hurstO contriDiversityO assymoriginalO]; %[time, only noise, rotation+noise, rotation+noise+transit, noise+transit]. 
%save('MultifIndex.txt','MultifIndex','-ascii'); 
f=figure;
surf(tauH,q,hurstO,'FaceAlpha',0.5,'FaceColor','interp');

colormap(f,jet);
colorbar
% Create xlabel
xlabel('\tau (days)','FontSize',12);

% Create ylabel
ylabel('q','FontSize',12);

% Create zlabel
zlabel('h(q,\tau)','FontSize',12);
zlim([0 1.5]);


