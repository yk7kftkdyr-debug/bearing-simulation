function [deltaw, clearInfo] = calcWorkingClearance(datafromvb)
%CALCWORKINGCLEARANCE Calculate roller-bearing working radial clearance.

Dw=datafromvb(2);Dm=datafromvb(3);deltar0=datafromvb(7);
e1=datafromvb(11);e2=datafromvb(12);
W1=datafromvb(17)*pi/30;W2=datafromvb(18)*pi/30;
o1=datafromvb(14);o2=datafromvb(15);
Di=datafromvb(35);Do=datafromvb(34);Ds=datafromvb(36);Dh=datafromvb(37);
es=datafromvb(38);eh=datafromvb(39);os=datafromvb(40);oh=datafromvb(41);
u1=datafromvb(42);u2=datafromvb(43);
Tr=datafromvb(44);To=datafromvb(45);Ti=datafromvb(46);Ta=datafromvb(47);
ruo1=datafromvb(48);ruo2=datafromvb(49);ruos=datafromvb(50);ruoh=datafromvb(51);
taos=datafromvb(52);taoh=datafromvb(53);tao1=datafromvb(54);tao2=datafromvb(55);taor=datafromvb(56);
Ts=datafromvb(57);Th=datafromvb(58);

Dr1=Dm+Dw+deltar0;
Dr2=Dm-Dw-deltar0;

if Ds==0
    deltai=2*u1*(Dr2/Di)/(((Dr2/Di)^2-1)*((((Dr2/Di)^2+1)/((Dr2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltai=2*u1*(Dr2/Di)/(((Dr2/Di)^2-1)*((((Dr2/Di)^2+1)/((Dr2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));
end
deltao=2*u2*(Do/Dr1)/(((Do/Dr1)^2-1)*((((Do/Dr1)^2+1)/((Do/Dr1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));

if deltai<=0
    deltai=0;
end
if deltao<=0
    deltao=0;
end
deltapd=-deltai-deltao;

deltat1=tao1*Do*(To-Ta);
deltat2=tao2*Di*(Ti-Ta);
deltatb=taor*Dw*(Tr-Ta);
deltats=taos*Di*(Ts-Ta);
deltath=taoh*Dh*(Th-Ta);
deltapt1=deltat1-2*deltatb-deltat2;
u1t=deltats-deltat2;
u2t=deltat1-deltath;

if deltai<=0 && (deltai+u1t)<=0
    u1t=0;
elseif deltai<=0 && (deltai+u1t)>0
    u1t=deltai+u1t;
elseif deltai>0 && (deltai+u1t)<=0
    u1t=-deltai;
end
if deltao<=0 && (deltao+u2t)<=0
    u2t=0;
elseif deltao<=0 && (deltao+u2t)>0
    u2t=deltao+u2t;
elseif deltao>0 && (deltao+u2t)<=0
    u2t=-deltao;
end

if Ds==0
    deltait=2*u1t*(Dr2/Di)/(((Dr2/Di)^2-1)*((((Dr2/Di)^2+1)/((Dr2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltait=2*u1t*(Dr2/Di)/(((Dr2/Di)^2-1)*((((Dr2/Di)^2+1)/((Dr2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));
end
deltaot=2*u2t*(Do/Dr1)/(((Do/Dr1)^2-1)*((((Do/Dr1)^2+1)/((Do/Dr1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));
deltapt2=-deltait-deltaot;
deltat=deltapt1+deltapt2;

Ri=(Dr2+Di)/4; Ro=(Dr1+Do)/4; Rs=(Ds+Di)/4; Rh=(Do+Dh)/4;
deltaf1=2*ruo2*Ri^3*W2^2/e2;
deltaf2=2*ruo1*Ro^3*W1^2/e1;
deltafs=2*ruos*Rs^3*W2^2/es;
deltafh=2*ruoh*Rh^3*W1^2/eh;
deltapf1=deltaf1-deltaf2;
u1f=deltafs-deltaf2;
u2f=deltaf1-deltafh;

if (deltai+deltats-deltat2)<=0 && (deltai+deltats-deltat2+u1f)<=0
    u1f=0;
elseif (deltai+deltats-deltat2)<=0 && (deltai+deltats-deltat2+u1f)>0
    u1f=deltai+deltats-deltat2+u1f;
elseif (deltai+deltats-deltat2)>0 && (deltai+deltats-deltat2+u1f)<=0
    u1f=-(deltai+deltats-deltat2);
end
if (deltao+deltat1-deltath)<=0 && (deltao+deltat1-deltath+u2f)<=0
    u2f=0;
elseif (deltao+deltat1-deltath)<=0 && (deltao+deltat1-deltath+u2f)>0
    u2f=deltao+deltat1-deltath+u2f;
elseif (deltao+deltat1-deltath)>0 && (deltao+deltat1-deltath+u2f)<=0
    u2f=-(deltao+deltat1-deltath);
end

if Ds==0
    deltaif=2*u1f*(Dr2/Di)/(((Dr2/Di)^2-1)*((((Dr2/Di)^2+1)/((Dr2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltaif=2*u1f*(Dr2/Di)/(((Dr2/Di)^2-1)*((((Dr2/Di)^2+1)/((Dr2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));
end
deltaof=2*u2f*(Do/Dr1)/(((Do/Dr1)^2-1)*((((Do/Dr1)^2+1)/((Do/Dr1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));
deltapf2=-deltaif-deltaof;
deltaf=deltapf1+deltapf2;
deltaw=deltar0+deltapd+deltat+deltaf;

clearInfo = struct('deltapd',deltapd,'deltapt1',deltapt1,'deltapt2',deltapt2, ...
    'deltapf1',deltapf1,'deltapf2',deltapf2,'deltat',deltat,'deltaf',deltaf,'deltaw',deltaw);
end
