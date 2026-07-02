function [returndata]=Famax1(datafromvb)         %计算球轴承接触区超出滚道的极限推力载荷

% 数据输入部分！！
n=datafromvb(1);Dw=datafromvb(2);Dm=datafromvb(3);f1=datafromvb(4);f2=datafromvb(5);a0=datafromvb(6)*pi/180;
Rp=datafromvb(7); %兜孔半径
ballden=datafromvb(61);e1=datafromvb(9);e2=datafromvb(10);e3=datafromvb(11);o1=datafromvb(12);o2=datafromvb(13);o3=datafromvb(14); 
W2=datafromvb(15)*pi/30;Fxx=datafromvb(16); Fyy=datafromvb(17); Fzz=datafromvb(18); Myy=datafromvb(19)+0.1; Mzz=datafromvb(20)+0.1; 
dangbianxishu1=datafromvb(21);dangbianxishu2=datafromvb(22);
cucao1=datafromvb(23);cucao2=datafromvb(24);cucao3=datafromvb(25); 
%润滑油基本参数
oilden=datafromvb(26); sita0=datafromvb(27);K=datafromvb(28); %常温下的导热系数
niandu0=datafromvb(29); nianya0=datafromvb(30); beita0=datafromvb(31);   %粘温系数（近似认为不变，wys论文中）   所有参数均为常温下的参数！！！
yindao=datafromvb(32); yindaojianxi=datafromvb(33);   % 引导方式,引导间隙
dianpianjiao=datafromvb(34);  %垫片角
%装配的相关参数
Do=datafromvb(36);Di=datafromvb(37);Ds=datafromvb(38);Dh=datafromvb(39);es=datafromvb(40);eh=datafromvb(41);os=datafromvb(42);oh=datafromvb(43);u1=datafromvb(44);u2=datafromvb(45);
taos=datafromvb(54);taoh=datafromvb(55);tao1=datafromvb(56);tao2=datafromvb(57);taob=datafromvb(58);Tb=datafromvb(46);To=datafromvb(47);Ti=datafromvb(48);Ta=datafromvb(49);Ts=datafromvb(59);Th=datafromvb(60);
ruo1=datafromvb(50);ruo2=datafromvb(51);ruos=datafromvb(52);ruoh=datafromvb(53);

W1=0; %外圈固定
%a0=acos(1-deltar0/(f1+f2-1)/Dw/2);
deltac0=Rp/2-Dw/2;
deltar0=2*(1-cos(a0))*(f1+f2-1)*Dw;               %求原始径向间隙！
deltaa0=2*sin(a0)*(f1+f2-1)*Dw ;                 %求原始轴向间隙！
Dy1=Dm+Dw+deltar0;Dy2=Dm-Dw-deltar0;  % 沟底直径
Dr1=Dy1-2*dangbianxishu1*Dw;Dr2=Dy2+2*dangbianxishu2*Dw;  % 挡边直径（引导面直径）  档边高系数给出！！！！

a0=acos(1-deltar0/(2*(f1+f2-1)*Dw)-(2*f2-1)*(1-cos(dianpianjiao*pi/180))/(2*(f1+f2-1)));  %存在垫片角ad时计算初始接触角

m=ballden*Dw*Dw*Dw*pi/6;  %球的质量计算
J=m*Dw*Dw/10;             %转动惯量计算
%求外、内圈赫兹接触载荷！
rou1=1/Dw*(4-1/f1-2*Dw*cos(a0)/(Dm+Dw*cos(a0)));
Rx1=Dw*(0.5*Dm/cos(a0)+0.5*Dw)/( Dm/cos(a0) );
Ry1=f1*Dw/(2*f1-1);
R21=Rx1*Ry1/(Rx1+Ry1);
E1=2/(((1-o1^2)/e1)+((1-o3^2)/e3));             %e1，e3为外圈和球的弹性模量，o1，o3为外圈，球的泊松比
K1=1.0339*(Ry1/Rx1)^0.636;
ee21=1.0003+0.5968/(Ry1/Rx1);          %ee21是外圈的ee2
ee11=((K1^2+1)*ee21-(K1^2-1)*ee21*( -2*cos(a0)/(Dm+Dw*cos(a0)) +1/(f1*Dw)  )/rou1)/2;

rou2=  1/Dw*(4-1/f2+2*Dw*cos(a0)/(Dm-Dw*cos(a0)));
E2=2/(((1-o2^2)/e2)+((1-o3^2)/e3));   %e1，e3为内圈和球的弹性模量，o1，o3为内圈，球的泊松比
Rx2=Dw*(0.5*Dm/cos(a0)-0.5*Dw)/(Dm/cos(a0) ); 
Ry2=f2*Dw/(2*f2-1);
K2=1.0339*(Ry2/Rx2)^0.636;
R22=Rx2*Ry2/(Rx2+Ry2);
R222=Dm/2+(f2-0.5)*Dw*cos(a0);
ee22=1.0003+0.5968/(Ry2/Rx2);  %ee22是内圈的ee2
ee12=((K2^2+1)*ee22-(K2^2-1)*ee22*( 2*cos(a0)/(Dm-Dw*cos(a0)) +1/(f2*Dw)  )/rou2)/2;

jiajiao2=acos(1-(Dr2-Dy2)/Dw);
a2=a0;  a21=a2+0.02;
K=16.7e8;  %查图6.5 Harris  仅与总曲率有关！！！

while (abs(a21-a2)>1e-6)
    a2=a21;
    a21=jiajiao2-asin( 0.0472*(2*K2^2*ee22/pi)^0.333*(cos(a0)/cos(a2)-1)^0.5/(Dw*rou2)^0.333);
end

Famax=n*Dw^2*K*sin(a2)*(cos(a0)/cos(a2)-1)^1.5;

%考虑配合以及热、离心力后的径向游隙 
if Ds==0
    deltai=2*u1*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltai=2*u1*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));  %Di,Do分别为轴承的内外圈直径;e2,es,o2,os分别是轴承内圈和轴的弹性模量以及泊松比;Ds表示空心轴的内径；u1表示轴与轴承内圈的过盈量
end
deltao=2*u2*(Do/Dy1)/(((Do/Dy1)^2-1)*((((Do/Dy1)^2+1)/((Do/Dy1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));  %eh,e1,oh,o1分别是轴承座和轴承外圈的弹性模量以及泊松比;Dh表示轴承座的外径；u2表示轴承座与轴承外圈的过盈量
if deltai<=0
    deltai0=deltai
    deltai=0
end
if deltao<=0
    deltao0=deltao
    deltao=0
end
deltapd=-deltai-deltao;  %装配中的过盈引起的游隙变化
%温度引起的膨胀
deltat1=tao1*Do*(To-Ta); %tao1,tao2,taob,taos,taoh表示轴承外圈、内圈、滚动体、轴以及轴承座材料的线性膨胀系数；To,Ti,Tb,Ts,Th,Ta分别表示轴承外圈、内圈、滚动体、轴、轴承座和环境的温度
deltat2=tao2*Di*(Ti-Ta);
deltatb=taob*Dw*(Tb-Ta);
deltats=taos*Di*(Ts-Ta);
deltath=taoh*Dh*(Th-Ta);
deltapt1=deltat1-2*deltatb-deltat2;%温度差引起的膨胀量
u1t=deltats-deltat2; u2t=deltat1-deltath;
if deltai<=0 && (deltai0+u1t)<=0
    u1t=0
elseif deltai<=0 && (deltai0+u1t)>0
    u1t=deltai0+u1t
elseif deltai>0 && (deltai+u1t)<=0
    u1t=-deltai
end
if deltao<=0 && (deltao0+u2t)<=0
    u2t=0
elseif deltao<=0 && (deltao0+u2t)>0
    u2t=deltao0+u2t
elseif deltao>0 && (deltao+u2t)<=0
    u2t=-deltao
end
if Ds==0
    deltait=2*u1t*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltait=2*u1t*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));  %Di,Do分别为轴承的内外圈直径;e2,es,o2,os分别是轴承内圈和轴的弹性模量以及泊松比;Ds表示空心轴的内径；u1表示轴与轴承内圈的过盈量
end
deltaot=2*u2t*(Do/Dy1)/(((Do/Dy1)^2-1)*((((Do/Dy1)^2+1)/((Do/Dy1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));  %eh,e1,oh,o1分别是轴承座和轴承外圈的弹性模量以及泊松比;Dh表示轴承座的外径；u2表示轴承座与轴承外圈的过盈量
deltapt2=-deltait-deltaot;  %温度差对装配量补偿引起的游隙变化
deltat=deltapt1+deltapt2;%温度引起的游隙变化
%离心力引起的变形
Ri=(Dy2+Dr2+2*Di)/8; %Dy,Dr代表沟道直径和挡边直径
Ro=(Dy1+Dr1+2*Do)/8;
Rs=(Ds+Di)/4;  Rh=(Do+Dh)/4;
deltaf1=2*ruo2*Ri^3*W2^2/e2;  %ruo1，ruo2代表外内圈材料密度
deltaf2=2*ruo1*Ro^3*W1^2/e1;
deltafs=2*ruos*Rs^3*W2^2/es;
deltafh=2*ruoh*Rh^3*W1^2/eh;
deltapf1=deltaf1-deltaf2;  %离心力引起的间隙
u1f=deltafs-deltaf2; u2f=deltaf1-deltafh;
if deltai<=0
    if (deltai0+deltats-deltat2)<=0 && (deltai0+deltats-deltat2+u1f)<=0
        u1f=0
    elseif (deltai0+deltats-deltat2)<=0 && (deltai0+deltats-deltat2+u1f)>0
        u1f=deltai0+deltats-deltat2+u1f
    elseif (deltai0+deltats-deltat2)>0 && (deltai0+deltats-deltat2+u1f)<=0
        u1f=-(deltai0+deltats-deltat2)
    end
else
    if (deltai+deltats-deltat2)<=0 && (deltai+deltats-deltat2+u1f)<=0
        u1f=0
    elseif (deltai+deltats-deltat2)<=0 && (deltai+deltats-deltat2+u1f)>0
        u1f=deltai+deltats-deltat2+u1f
    elseif (deltai+deltats-deltat2)>0 && (deltai+deltats-deltat2+u1f)<=0
        u1f=-(deltai+deltats-deltat2)
    end
end
if deltao<=0
    if (deltao0+deltat1-deltath)<=0 && (deltao0+deltat1-deltath+u2f)<=0
        u2f=0
    elseif (deltao0+deltat1-deltath)<=0 && (deltao0+deltat1-deltath+u2f)>0
        u2f=deltao0+deltat1-deltath+u2f
    elseif (deltao0+deltat1-deltath)>0 && (deltao0+deltat1-deltath+u2f)<=0
        u2f=-(deltao0+deltat1-deltath)
    end
else
    if (deltao+deltat1-deltath)<=0 && (deltao+deltat1-deltath+u2f)<=0
        u2f=0
    elseif (deltao+deltat1-deltath)<=0 && (deltao+deltat1-deltath+u2f)>0
        u2f=deltao+deltat1-deltath+u2f
    elseif (deltao+deltat1-deltath)>0 && (deltao+deltat1-deltath+u2f)<=0
        u2f=-(deltao+deltat1-deltath)
    end
end
if Ds==0
    deltaif=2*u1f*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(1-os)));
else
    deltaif=2*u1f*(Dy2/Di)/(((Dy2/Di)^2-1)*((((Dy2/Di)^2+1)/((Dy2/Di)^2-1)+o2)+e2/es*(((Di/Ds)^2+1)/((Di/Ds)^2-1)-os)));  %Di,Do分别为轴承的内外圈直径;e2,es,o2,os分别是轴承内圈和轴的弹性模量以及泊松比;Ds表示空心轴的内径；u1表示轴与轴承内圈的过盈量
end
deltaof=2*u2f*(Do/Dy1)/(((Do/Dy1)^2-1)*((((Do/Dy1)^2+1)/((Do/Dy1)^2-1)-o1)+e1/eh*(((Dh/Do)^2+1)/((Dh/Do)^2-1)+oh)));  %eh,e1,oh,o1分别是轴承座和轴承外圈的弹性模量以及泊松比;Dh表示轴承座的外径；u2表示轴承座与轴承外圈的过盈量
deltapf2=-deltaif-deltaof;  %离心力对装配量补偿引起的游隙变化
deltaf=deltapf1+deltapf2;%离心力引起的游隙变化
deltaw=deltar0+deltapd+deltat+deltaf;%工作游隙
save Famax;
save deltapd;save deltapt1;save deltapt2;save deltapf1;save deltapf2;save deltaw



