function [deltaw, clearInfo] = calcBallWorkingClearance(datafromvb)
%CALCBALLWORKINGCLEARANCE Calculate ball-bearing working radial clearance.

Dw=datafromvb(2);Dm=datafromvb(3);f1=datafromvb(4);f2=datafromvb(5);a0=datafromvb(6)*pi/180;
deltar0=2*(1-cos(a0))*(f1+f2-1)*Dw;
Dy1=Dm+Dw+deltar0;Dy2=Dm-Dw-deltar0;

e1=datafromvb(9);e2=datafromvb(10);
o1=datafromvb(12);o2=datafromvb(13);
W1=datafromvb(69)*pi/30;W2=datafromvb(15)*pi/30;
Do=datafromvb(36);Di=datafromvb(37);Ds=datafromvb(38);Dh=datafromvb(39);
es=datafromvb(40);eh=datafromvb(41);os=datafromvb(42);oh=datafromvb(43);
u1=datafromvb(44);u2=datafromvb(45);
Tb=datafromvb(46);To=datafromvb(47);Ti=datafromvb(48);Ta=datafromvb(49);
ruo1=datafromvb(50);ruo2=datafromvb(51);ruos=datafromvb(52);ruoh=datafromvb(53);
taos=datafromvb(54);taoh=datafromvb(55);tao1=datafromvb(56);tao2=datafromvb(57);taob=datafromvb(58);
Ts=datafromvb(59);Th=datafromvb(60);

if Ds==0
    deltai=2*u1*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltai=2*u1*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));
end
deltao=2*u2*(Do/Dy1)/(((Do/Dy1)^2-1)*((((Do/Dy1)^2+1)/((Do/Dy1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));

deltai0=deltai;
deltao0=deltao;
if deltai<=0
    deltai=0;
end
if deltao<=0
    deltao=0;
end
deltapd=-deltai-deltao;

deltat1=tao1*Do*(To-Ta);
deltat2=tao2*Di*(Ti-Ta);
deltatb=taob*Dw*(Tb-Ta);
deltats=taos*Di*(Ts-Ta);
deltath=taoh*Dh*(Th-Ta);
deltapt1=deltat1-2*deltatb-deltat2;
u1t=deltats-deltat2;
u2t=deltat1-deltath;

if deltai<=0 && (deltai0+u1t)<=0
    u1t=0;
elseif deltai<=0 && (deltai0+u1t)>0
    u1t=deltai0+u1t;
elseif deltai>0 && (deltai+u1t)<=0
    u1t=-deltai;
end
if deltao<=0 && (deltao0+u2t)<=0
    u2t=0;
elseif deltao<=0 && (deltao0+u2t)>0
    u2t=deltao0+u2t;
elseif deltao>0 && (deltao+u2t)<=0
    u2t=-deltao;
end

if Ds==0
    deltait=2*u1t*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltait=2*u1t*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));
end
deltaot=2*u2t*(Do/Dy1)/(((Do/Dy1)^2-1)*((((Do/Dy1)^2+1)/((Do/Dy1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));
deltapt2=-deltait-deltaot;
deltat=deltapt1+deltapt2;

Ri=(Dy2+Di)/4; Ro=(Dy1+Do)/4; Rs=(Ds+Di)/4; Rh=(Do+Dh)/4;
deltaf1=2*ruo2*Ri^3*W2^2/e2;
deltaf2=2*ruo1*Ro^3*W1^2/e1;
deltafs=2*ruos*Rs^3*W2^2/es;
deltafh=2*ruoh*Rh^3*W1^2/eh;
deltapf1=deltaf1-deltaf2;
u1f=deltafs-deltaf2;
u2f=deltaf1-deltafh;

if deltai<=0
    fiti=deltai0+deltats-deltat2;
else
    fiti=deltai+deltats-deltat2;
end
if deltao<=0
    fito=deltao0+deltat1-deltath;
else
    fito=deltao+deltat1-deltath;
end

if fiti<=0 && (fiti+u1f)<=0
    u1f=0;
elseif fiti<=0 && (fiti+u1f)>0
    u1f=fiti+u1f;
elseif fiti>0 && (fiti+u1f)<=0
    u1f=-fiti;
end
if fito<=0 && (fito+u2f)<=0
    u2f=0;
elseif fito<=0 && (fito+u2f)>0
    u2f=fito+u2f;
elseif fito>0 && (fito+u2f)<=0
    u2f=-fito;
end

if Ds==0
    deltaif=2*u1f*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltaif=2*u1f*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));
end
deltaof=2*u2f*(Do/Dy1)/(((Do/Dy1)^2-1)*((((Do/Dy1)^2+1)/((Do/Dy1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));
deltapf2=-deltaif-deltaof;
deltaf=deltapf1+deltapf2;
deltaw=deltar0+deltapd+deltat+deltaf;

clearInfo = struct('deltapd',deltapd,'deltapt1',deltapt1,'deltapt2',deltapt2, ...
    'deltapf1',deltapf1,'deltapf2',deltapf2,'deltat',deltat,'deltaf',deltaf,'deltaw',deltaw);
end
