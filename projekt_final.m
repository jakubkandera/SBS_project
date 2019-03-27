%% Mikrovlnný smìrový spoj
close all
clear all 

%% 1.skok

f=10*10^9;
Gant=40;%1.20 m gain mid band (dBi
EbN = 10.5;%QPSK
% EbN = 15;%16QAM
fb=40*10^6; %56Mb/s
B=28*10^6;%40MHz


c=3*10^8;

d1=4.79*10^3;
EIRP=20;%dBm

FSL=((4*pi*d1*f)/c)^2;
FSLdB=10*log10(FSL);
Lgas=0.05*4.79;
Lcable=0.6;
Pp1=EIRP-FSLdB+Gant-Lcable-Lgas;

Re=6371*10^3;
k1=4/3;
k2=0.78;
Reff=k1*Re;
Reff2=k2*Re;
d01=2.63*10^3;
d02=d1-2.63*10^3;
x=(d01*d02)/(2*Reff);
x2=(d01*d02)/(2*Reff2);


F1=17.3* sqrt(((d01*d02)*10^-6)/(f*10^-9*d1*10^-3));

%% 2.skok

d2=10.59*10^3;
FSL=((4*pi*d2*f)/c)^2;
FSLdB=10*log10(FSL);
Lgas=0.05*14;
Lcable=0.6;
Pp2=EIRP-FSLdB+Gant-Lcable-Lgas;
% 
% Re=6371*10^3;
% k1=4/3;
% k2=0.78;
% Reff=k1*Re;
% Reff2=k2*Re;
% d01=5.3*10^3;
% d02=d-5.3*10^3;
% x3=(d1*d2)/(2*Reff);
% x4=(d1*d2)/(2*Reff2);
% h2=429;
% 
% F2=17.3* sqrt(((d01*d02)*10^-6)/(f*10^-9*d*10^-3));

%% 3skok

d3=11.09*10^3;

FSL=((4*pi*d3*f)/c)^2;
FSLdB=10*log10(FSL);
Lgas=0.05*7.36;
Lcable=0.6;
Pp3=EIRP-FSLdB+Gant-Lcable-Lgas;

% Re=6371*10^3;
% k1=4/3;
% k2=0.78;
% Reff=k1*Re;
% Reff2=k2*Re;
% d01=1.5*10^3;
% d02=d-1.5*10^3;
% x5=(d01*d02)/(2*Reff);
% x6=(d01*d02)/(2*Reff2);
% 
% F3=17.3* sqrt(((d01*d02)*10^-6)/(f*10^-9*d*10^-3));


%% 4sko

d4=2.87*10^3;

FSL=((4*pi*d4*f)/c)^2;
FSLdB=10*log10(FSL);
Lgas=0.05*2.87;
Lcable=0.6;
Pp4=EIRP-FSLdB+Gant-Lcable-Lgas;

Re=6371*10^3;
k1=4/3;
k2=0.78;
Reff=k1*Re;
Reff2=k2*Re;
d01=1.5*10^3;
d02=d1-1.5*10^3;
x7=(d01*d02)/(2*Reff);
x8=(d01*d02)/(2*Reff2);

F4=17.3* sqrt(((d01*d02)*10^-6)/(f*10^-9*d1*10^-3));


%% %citlivost pøijímaèe
% EbN = 10.5;%QPSK
% EbN = 15;%16QAM
% fb=17*10^6; %56Mb/s
% B=13.75*10^6;%40MHz
% B=30.5*10^6;%40MHz
k=1.38*10^-23;
T=7+273.15;

N=k*T*B;
NdB=10*log10((1*10^3)*N);
%CN=EbN*(fb/B);
%CNdB=10*log10(CN);
CNdB=EbN + 10*log10(fb/B);
%SOR=CN*N;
%SORdB=10*log10(SOR*1000);%dBm
SORdB=CNdB + NdB;%dBm


%% Margin
margin1=Pp1-SORdB;
margin2=Pp2-SORdB;
margin3=Pp3-SORdB;
margin4=Pp4-SORdB;

%% Fading and enhancement due to multipath and related mechanisms

%pro skok1
d1=4.79;
dN1=-200;
K=1*10^(-4.6-0.0027*dN1);
hr=30+451;%vyska vysilaèe + nadmoøska vyska
he=14+275;

ep=(abs(hr-he))/d1; % hr- vý¹ka pøijímaèe, he- vý¹ka vysílaèe
epabs=abs(ep);
hl=min(he,hr);
A1=margin1;
pw=K*(d1^3.1) * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl-A1/10); %A=margin
Pns1=pw/100;

%pro skok 2
d2=10.59;
% dN1=-200;
K=1*10^(-4.6-0.0027*dN1);
hr=30+400;
he=30+451;

ep=(abs(hr-he))/d2;
epabs=abs(ep);
hl=min(he,hr);
A2=margin2;
pw=K*d2^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl-A2/10);
Pns2=pw/100;

%% pro skok 3
d3=11.09;
% dN1=-200;
K=1*10^(-4.6-0.0027*dN1);
hr=63+356;
he=30+400;

ep=(abs(hr-he))/d3;
epabs=abs(ep);
hl=min(he,hr);
A3=margin3;
pw=K*d3^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl-A3/10);
Pns3=pw/100;

%pro skok4
d4=2.87;
% dN1=-400;% hodnota mezi -200 a¾ -400
K=1*10^(-4.6-0.0027*dN1);
hr=40+215;
he=63+356;

ep=(abs(hr-he))/d4;
epabs=abs(ep);
hl=min(he,hr);
A4=margin4;
pw=K*d4^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl-A4/10);


Pns4=pw/100;


Pns=[Pns1 Pns2 Pns3 Pns4];

Psum=0;
for i=1:4
    
Psum=Psum+Pns(i);

end

A=[A1 A2 A3 A4];
c=0;
dsum=0;
P=0;
P1=0;
d=[d1 d2 d3 d4];
Psuma=0;
for n=1:3
%dsum=dsum+d(n);
c(n)=0.5+0.0052*A(n)+0.0025*(d(n)+d(n+1));
%P1=P1+Pns(n);
%P(n)=(Pns(n)*P1)^c(n);
Psuma=Psuma+((Pns(n)*Pns(n+1))^c(n));

end

Pvice=Psum-Psuma


%% utlum destem

R=32; %mereni CHMU
%pro 10Ghz

kv=0.01129;% horizontalni polarizace nema takovy vliv, muzu zanedbat
av=1.2156;
kapa=kv/2;
alfa=(kv*av)/(2*kapa);
gama= kapa*R^alfa;

def1=d1/(1+(d1/(35*exp(-0.015*R))));
def2=d2/(1+(d2/(35*exp(-0.015*R))));
def3=d3/(1+(d3/(35*exp(-0.015*R))));
def4=d4/(1+(d4/(35*exp(-0.015*R))));

r1=1/(0.477*d1^0.633 *R^(0.073*alfa)*(f*10^-9)^0.123 - 10.579*(1-exp(-0.024*d1)));
r2=1/(0.477*d2^0.633 *R^(0.073*alfa)*(f*10^-9)^0.123 - 10.579*(1-exp(-0.024*d2)));
r3=1/(0.477*d3^0.633 *R^(0.073*alfa)*(f*10^-9)^0.123 - 10.579*(1-exp(-0.024*d3)));
r4=1/(0.477*d4^0.633 *R^(0.073*alfa)*(f*10^-9)^0.123 - 10.579*(1-exp(-0.024*d4)));

A01=(gama*d1*r1);
A02=(gama*d2*r2);
A03=(gama*d3*r3);
A04=(gama*d4*r4);

if (f*10^-9)>=10
    
    c0=0.12+0.4*(log10(f*10^-9/10)^0.8);
else
  
    c0=0.12;
end
c1=(0.07^c0)*(0.12^(1-c0));
c2=0.855*c0+0.546*(1-c0);
c3=0.139*c0+0.043*(1-c0);
p=[0:0.00001:1];

Ap1=(A01*c1*p.^-(c2+c3.*log10(p)));
Ap2=(A02*c1*p.^-(c2+c3.*log10(p)));
Ap3=(A03*c1*p.^-(c2+c3.*log10(p)));
Ap4=(A04*c1*p.^-(c2+c3.*log10(p)));

Pdest1=0; % rezerva unik je dostatecne velika (25dB)
Pdest2=0;
Pdest3=0;
Pdest4=0; % rezerva unik je dostatecne velika (30dB)

Pdest=Pdest2

figure(1);
plot(Ap1, p);
set(gca, 'YScale', 'log')
grid on;
title('Ztraty vlivem deste- pro skok 1 ');
ylabel('p[-]');
xlabel('A[dB]');

figure(2);
plot(Ap2, p);
set(gca, 'YScale', 'log')
grid on;
title('Ztraty vlivem deste- pro skok 2 ');
ylabel('p[-]');
xlabel('A[dB]');

figure(3);
plot(Ap3, p);
set(gca, 'YScale', 'log')
grid on;
title('Ztraty vlivem deste- pro skok 3 ');
ylabel('p[-]');
xlabel('A[dB]');


figure(4);
plot(Ap4, p);
set(gca, 'YScale', 'log')
grid on;
title('Ztraty vlivem deste- pro skok 4 ');
ylabel('p[-]');
xlabel('A[dB]');


%% Reduction of cross-polar discrimination (XPD)
pw01=K*d1^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl); %A=margin
P01=pw01/100;
pw02=K*d2^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl); %A=margin
P02=pw02/100;
pw03=K*d3^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl); %A=margin
P03=pw03/100;
pw04=K*d4^3.1 * (1+epabs)^-1.29 *(f*10^-9)^0.8 * 10^(-0.00089*hl); %A=margin
P04=pw04/100;



ucinnost01= 1-exp(-0.2*(P01)^0.75);
ucinnost03= 1-exp(-0.2*(P03)^0.75);
ucinnost02= 1-exp(-0.2*(P02)^0.75);
ucinnost04= 1-exp(-0.2*(P04)^0.75);
Q01=-10*log10((0.7*ucinnost01)/P01);
Q02=-10*log10((0.7*ucinnost02)/P02);
Q03=-10*log10((0.7*ucinnost03)/P03);
Q04=-10*log10((0.7*ucinnost04)/P04);

XPD0=28+5;%nalezeno ve specifikaci anteny

C01=XPD0+Q01;
C02=XPD0+Q02;
C03=XPD0+Q03;
C04=XPD0+Q04;

Mxpd01=C01-CNdB;
Mxpd02=C02-CNdB;
Mxpd03=C03-CNdB;
Mxpd04=C04-CNdB;

Pxp01=P01*10^(-Mxpd01/10);
Pxp02=P02*10^(-Mxpd02/10);
Pxp021=P03*10^(-Mxpd03/10);
Pxp04=P04*10^(-Mxpd04/10);

Pxp=Pxp01+Pxp02+Pxp021+Pxp04
%% Distortion due to propagation effects (inter symbol interference, group delay)
Kn=1; %jestli teda mame QPSK
tm1=(10^-9)*0.7*(d1/50)^1.3 ;
tm2=(10^-9)*0.7*(d2/50)^1.3 ;
tm3=(10^-9)*0.7*(d3/50)^1.3 ;
tm4=(10^-9)*0.7*(d4/50)^1.3 ;
m=2;
SR=fb/(m*(3/4));
T=1/SR;
Ps1= 2.15*ucinnost01*Kn*(tm1/T)^2;
Ps2= 2.15*ucinnost02*Kn*(tm2/T)^2;
Ps3= 2.15*ucinnost03*Kn*(tm3/T)^2;
Ps4= 2.15*ucinnost04*Kn*(tm4/T)^2;


Ps=Ps1+Ps2+Ps3+Ps4


%% Total outage prediction

pt=Pvice+Pxp+Ps+Pdest
ptperc=pt*100



