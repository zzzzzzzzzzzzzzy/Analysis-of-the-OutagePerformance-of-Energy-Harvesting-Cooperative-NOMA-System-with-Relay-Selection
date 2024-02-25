clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOMA分配系数变化图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = 0.1;
num = 20;
R0=2;
T=2^R0-1;
xx=T;
N=1;

hsr = 1;
hrd = 1;
hsd = 1;
Ps=10.^((0:30)./10);

a2 = 1-a1;

[x,y] = bosong(num);
dsr = sqrt(x.^2+y.^2)./sqrt(2);
drd = sqrt((1-x).^2+(1-y).^2)./sqrt(2);
dsd = 1;


alpha = 0.1;
beta = 0.1;

for i = 1:31
    P = Ps(i);
    G = sqrt(P./(P.*hsr.^2./(1.+dsr)+1));
    gamma_sr = hsr.^2.*P./(1+dsr);
    gamma_sd = P.*hsd.^2.*a1./(P.*a2.*hsd.^2+1+dsd);
    gamma_rd = (P.*G.^2.*hrd.^2.*hsr.^2.*P.*a1.*(1-alpha-beta))./(P.*G.^2.*hrd.^2.*hsr.^2.*P.*a2.*(1-alpha-beta)+G.^2.*hrd.^2.*(1+dsr)+(1-drd).*(1+dsr));

    for j = 1:length(gamma_rd)-1
        a(j) = min(gamma_rd(j),min(gamma_sd,gamma_sr(j)));
    end
    [gama,index2] = max(a);
    
    [gamma,index1] = max(gamma_rd.*gamma_sr./(gamma_sr+gamma_rd+1));

    SD = 1./gamma_sd;
    RR = 1./0;
    syms yy;

    SR1 = 1./gamma_sr(index1);
    RD1 = 1./gamma_rd(index1);
    SR2 = 1./gamma_sr(index2);
    RD2 = 1./gamma_rd(index2);

    z1(i) = (1-(RD1.*exp(-RD1.*xx)).*(int((exp((-((SR1.*(xx+yy+1).*xx)./yy))-(RD1.*yy)))./(1+((SR1.*(xx+yy+1)*xx)./(RR.*yy))),yy,0,inf))).^N;
    z2(i) = (1-(RD2.*exp(-RD2.*xx)).*(int((exp((-((SR2.*(xx+yy+1).*xx)./yy))-(RD2.*yy)))./(1+((SR2.*(xx+yy+1)*xx)./(RR.*yy))),yy,0,inf))).^N;

end
semilogy(10*log10(Ps),z1,'--gS')
grid on 
hold on
semilogy(10*log10(Ps),z2,'S-g')
grid on 
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = 0.55;
num = 20;
R0=2;
T=2^R0-1;
xx=T;
N=1;

hsr = 1;
hrd = 1;
hsd = 1;
Ps=10.^((0:30)./10);

a2 = 1-a1;

[x,y] = bosong(num);
dsr = sqrt(x.^2+y.^2)./sqrt(2);
drd = sqrt((1-x).^2+(1-y).^2)./sqrt(2);
dsd = 1;


alpha = 0.1;
beta = 0.1;

for i = 1:31
    P = Ps(i);
    G = sqrt(P./(P.*hsr.^2./(1.+dsr)+1));
    gamma_sr = hsr.^2.*P./(1+dsr);
    gamma_sd = P.*hsd.^2.*a1./(P.*a2.*hsd.^2+1+dsd);
    gamma_rd = (P.*G.^2.*hrd.^2.*hsr.^2.*P.*a1.*(1-alpha-beta))./(P.*G.^2.*hrd.^2.*hsr.^2.*P.*a2.*(1-alpha-beta)+G.^2.*hrd.^2.*(1+dsr)+(1-drd).*(1+dsr));

    for j = 1:length(gamma_rd)-1
        a(j) = min(gamma_rd(j),min(gamma_sd,gamma_sr(j)));
    end
    [gama,index2] = max(a);
    
    [gamma,index1] = max(gamma_rd.*gamma_sr./(gamma_sr+gamma_rd+1));

    SD = 1./gamma_sd;
    RR = 1./0;
    syms yy;

    SR1 = 1./gamma_sr(index1);
    RD1 = 1./gamma_rd(index1);
    SR2 = 1./gamma_sr(index2);
    RD2 = 1./gamma_rd(index2);

    z1(i) = (1-(RD1.*exp(-RD1.*xx)).*(int((exp((-((SR1.*(xx+yy+1).*xx)./yy))-(RD1.*yy)))./(1+((SR1.*(xx+yy+1)*xx)./(RR.*yy))),yy,0,inf))).^N;
    z2(i) = (1-(RD2.*exp(-RD2.*xx)).*(int((exp((-((SR2.*(xx+yy+1).*xx)./yy))-(RD2.*yy)))./(1+((SR2.*(xx+yy+1)*xx)./(RR.*yy))),yy,0,inf))).^N;

end
semilogy(10*log10(Ps),z1,'--bo')
grid on 
hold on
semilogy(10*log10(Ps),z2,'bo-')
grid on 
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = 0.6;
num = 20;
R0=2;
T=2^R0-1;
xx=T;
N=1;

hsr = 1;
hrd = 1;
hsd = 1;
Ps=10.^((0:30)./10);

a2 = 1-a1;

[x,y] = bosong(num);
dsr = sqrt(x.^2+y.^2)./sqrt(2);
drd = sqrt((1-x).^2+(1-y).^2)./sqrt(2);
dsd = 1;


alpha = 0.1;
beta = 0.1;

for i = 1:31
    P = Ps(i);
    G = sqrt(P./(P.*hsr.^2./(1.+dsr)+1));
    gamma_sr = hsr.^2.*P./(1+dsr);
    gamma_sd = P.*hsd.^2.*a1./(P.*a2.*hsd.^2+1+dsd);
    gamma_rd = (P.*G.^2.*hrd.^2.*hsr.^2.*P.*a1.*(1-alpha-beta))./(P.*G.^2.*hrd.^2.*hsr.^2.*P.*a2.*(1-alpha-beta)+G.^2.*hrd.^2.*(1+dsr)+(1-drd).*(1+dsr));

    for j = 1:length(gamma_rd)-1
        a(j) = min(gamma_rd(j),min(gamma_sd,gamma_sr(j)));
    end
    [gama,index2] = max(a);
    
    [gamma,index1] = max(gamma_rd.*gamma_sr./(gamma_sr+gamma_rd+1));

    SD = 1./gamma_sd;
    RR = 1./0;
    syms yy;

    SR1 = 1./gamma_sr(index1);
    RD1 = 1./gamma_rd(index1);
    SR2 = 1./gamma_sr(index2);
    RD2 = 1./gamma_rd(index2);
    z1(i) = (1-(RD1.*exp(-RD1.*xx)).*(int((exp((-((SR1.*(xx+yy+1).*xx)./yy))-(RD1.*yy)))./(1+((SR1.*(xx+yy+1)*xx)./(RR.*yy))),yy,0,inf))).^N;
    z2(i) = (1-(RD2.*exp(-RD2.*xx)).*(int((exp((-((SR2.*(xx+yy+1).*xx)./yy))-(RD2.*yy)))./(1+((SR2.*(xx+yy+1)*xx)./(RR.*yy))),yy,0,inf))).^N;

end
semilogy(10*log10(Ps),z1,'--m*')
grid on 
hold on
semilogy(10*log10(Ps),z2,'m*-')
grid on 
hold on
correct([0.865,0.8875,0.91,0.9325,0.955,0.9775],{'10^-^3','10^-^2^.^5','10^-^2','10^-^1^.^5','10^-^1','10^-^0^.^5'});


legend('F1 a=0.5','F2 a=0.5','F1 a=0.55','F2 a=0.55','F1 a=0.6','F2 a=0.6');
xlabel('Transmit Power [dB]') 
ylabel('Outage Probability')

