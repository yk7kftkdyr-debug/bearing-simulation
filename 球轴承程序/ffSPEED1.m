function ffSPEED1(datafromvb,ttspeed)   % 求Wo Wx Wz

load loadi;load www3; loadj=length(loadi);
ballspeed=ttspeed;

for i=1:loadj
    x(i)=www3(i);   r(i)=www3(loadj+i);
end 
ffLOAD(datafromvb,loadi,www3); load q1q2a1a2;load Ph1;load Ph2;load aa1;load aa2
for i=1:loadj
    Q1(i)=q1q2a1a2(i);    Q2(i)=q1q2a1a2(loadj+i);  a1(i)=q1q2a1a2(2*loadj+i);  a2(i)=q1q2a1a2(3*loadj+i);
end
for i=1:loadj
    Wo(i)=ballspeed(i);      Wx(i)=ballspeed(loadj+i);    Wz(i)=ballspeed(2*loadj+i);  Wy(i)=-Wo(i)/100;%计算所得
end

if loadj==datafromvb(1)
    Wononload=0; Wxnonload=0;
else
Wononload=ballspeed(3*loadj+1); Wxnonload=ballspeed(3*loadj+2); Wznonload=0; Wynonload=0;%初始值
end

for i=1:loadj
    if Wz(i)<0
        if Wz(n+1-i)>0
            Wz(i)=Wz(n+1-i);
        else
            Wz(i)=abs(Wz(i));
        end
    end  
end

Wc=(sum(Wo)+ Wononload*(n-loadj))/n;         % 保持架转速为各滚动体公转速度的平均值！！


for i=1:loadj
    if Q1(i)==max(Q1)
       markQ1=i;          % 找到受载最大的球
    end
end


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
lubricationtype=datafromvb(70); %轴承润滑形式
fcbr=datafromvb(71); fcbc=datafromvb(72); fccr=datafromvb(73);  %球/滚道、球/保持架、保持架/套圈之间固体润滑摩擦系数

deltac0=Rp/2-Dw/2;
deltar0=2*(1-cos(a0))*(f1+f2-1)*Dw;               %求原始径向间隙！
deltaa0=2*sin(a0)*(f1+f2-1)*Dw ;                 %求原始轴向间隙！
Dy1=Dm+Dw+deltar0;Dy2=Dm-Dw-deltar0;  % 沟底直径
Dr1=Dy1-2*dangbianxishu1*Dw;Dr2=Dy2+2*dangbianxishu2*Dw;  % 挡边直径（引导面直径）  档边高系数给出！！！！

a0=acos(1-deltar0/(2*(f1+f2-1)*Dw)-(2*f2-1)*(1-cos(dianpianjiao*pi/180))/(2*(f1+f2-1)));  %存在垫片角ad时计算初始接触角

z7=0;z8=0;z9=0;z10=0;z11=0;z12=0;z13=0;z14=0;z15=0; i=1;j=1;delt=0;
      
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
ee22=1.0003+0.5968/(Ry2/Rx2);  %ee21是内圈的ee2
ee12=((K2^2+1)*ee22-(K2^2-1)*ee22*( 2*cos(a0)/(Dm-Dw*cos(a0)) +1/(f2*Dw)  )/rou2)/2;


for i=1:loadj
sita(i)=2*pi*(loadi(i))/n;       % 从1号球开始考虑，与滚子轴承程序不同！！
if lubricationtype==0
    % 黏度计算！！
    %relands粘压关系，王艳霜博士论文P35
    %ZZ=nianya0/(5.1*1e-9)/(log(niandu0)+9.67);
    %这里不考虑粘温关系如下：
    %S0=beita0*(sita0+273-138)/(log(niandu0)+9.67);
    
    
    oilden1(i)=oilden*(1+(0.6e-9*Ph1(i))/(1+1.7e-9*Ph1(i)) );   % 滑油密度计算    公式在清华弹流12页 
    niandu01(i)=0.132*sinh(925.19/(-148.85+sita0+273))*oilden1(i)/1e6;          %sita0度时的粘度，wys第40页
    niandu1(i)=niandu01(i)*exp(Ph1(i)*nianya0);      % 考虑粘压系数后的 黏度！ 
    
    oilden2(i)=oilden*(1+(0.6e-9*Ph2(i))/(1+1.7e-9*Ph2(i)) );   % 滑油密度计算    公式在清华弹流12页 
    niandu02(i)=0.132*sinh(925.19/(-148.85+sita0+273))*oilden2(i)/1e6; %120度时的粘度，wys第40页
    niandu2(i)=niandu02(i)*exp(Ph2(i)*nianya0);      % 考虑粘压系数后的 黏度！
    
    %niandu2(i)=niandu0*exp((log(niandu0)+9.67)*(-1+(5.1*1e-9*Ph2(i))*exp(ZZ)));
    
    %求拖动系数%求球沿椭圆短轴轴方向的摩擦力！
    Re=2*f1*Dw/(2*f1+1);   
    deltaU11(i)=abs(-Wo(i)*Dm/2+(abs(Wx(i))*cos(a1(i))+ abs(Wz(i))*sin(a1(i))-Wo(i)*cos(a1(i)))*(Re-(Re^2-aa1(i)^2)^0.5+((Dw/2)^2-aa1(i)^2)^0.5));
    U11(i)=abs(-Wo(i)*Dm/4-0.5*(abs(Wx(i))*cos(a1(i))+ abs(Wz(i))*sin(a1(i))+Wo(i)*cos(a1(i)))*(Re-(Re^2-aa1(i)^2)^0.5+((Dw/2)^2-aa1(i)^2)^0.5));
    S11(i)=0;%deltaU11(i)/U11(i);
    deltaU1(i)=abs(-Wo(i)*Dm/2+(abs(Wx(i))*cos(a1(i))+ abs(Wz(i))*sin(a1(i))-Wo(i)*cos(a1(i)))*((Re^2-(aa1(i))^2)^0.5-(Re^2-aa1(i)^2)^0.5+((Dw/2)^2-aa1(i)^2)^0.5));
    U1(i)=abs(-Wo(i)*Dm/4-0.5*(abs(Wx(i))*cos(a1(i))+ abs(Wz(i))*sin(a1(i))+Wo(i)*cos(a1(i)))*((Re^2-(aa1(i))^2)^0.5-(Re^2-aa1(i)^2)^0.5+((Dw/2)^2-aa1(i)^2)^0.5));
    S12(i)=abs(deltaU1(i)/U1(i));
    S1(i)=(S11(i)+S12(i))/2;
    P01(i)=Q1(i)/Rx1^2/E1;
    U01(i)=niandu1(i)*U1(i)/(E1*Rx1);
    sita1(i)=sita0*(K/niandu1(i)/(U1(i)^2));
    A1(i)=-1.3527*10^(-20)*P01(i)^1.0849*U01(i)^(-0.28538)*sita1(i)^(-1.8235);
    B1(i)=1.0217*10^(-25)*P01(i)^2.8202*U01(i)^(-0.53722)*sita1(i)^(-3.7184);
    C1(i)=1.6836*10^(-23)*P01(i)^1.7563*U01(i)^(-0.25749)*sita1(i)^(-2.6437);
    D1(i)=8.0232*10^(-23)*P01(i)^1.1072*U01(i)^(-0.3333)*sita1(i)^(-2.0078);
    
   % A1(i)=-0.0196;
   % B1(i)=0.39873;
   % C1(i)=39.002;
   % D1(i)=0.01984;
    
    miu1(i)=0.0127*(50/(50-S1(i)))*log(0.584*Q1(i)/niandu0/deltaU1(i)/(U1(i))^2);  % 滚动体/内圈的拖动系数!%CHENZHE GAI Bennit公式
    T1(i)=miu1(i)*Q1(i);                      % 滚动体/外圈的拖动力！
    
    Ri=2*f2*Dw/(2*f2+1); % Ri=Rx2; 
    %y方向相对滑动
    deltaU21(i)=abs(Wy(i)*(Ri-(Ri^2-aa2(i)^2)^0.5+((Dw/2)^2-aa2(i)^2)^0.5));
    U21(i)=0.5*abs(Wy(i)*(Ri-(Ri^2-aa2(i)^2)^0.5+((Dw/2)^2-aa2(i)^2)^0.5));
    S21(i)=0;
   %滑动速度deltaU2=Ub-Ur
    deltaU2(i)=abs(+(W2-Wo(i))*Dm/2-( abs(Wx(i))*cos(a2(i))+ abs (Wz(i))*sin(a2(i))+(W2-Wo(i))*cos(a2(i)))*((Ri^2-(aa2(i))^2)^0.5-(Ri^2-aa2(i)^2)^0.5+((Dw/2)^2-aa2(i)^2)^0.5));
   %卷吸速度U2=(Ub+Ur)/2
    U2(i)=abs(+(W2-Wo(i))*Dm/4+0.5*( abs(Wx(i))*cos(a2(i))+abs(Wz(i))*sin(a2(i))-(W2-Wo(i))*cos(a2(i)))*((Ri^2-(aa2(i))^2)^0.5-(Ri^2-aa2(i)^2)^0.5+((Dw/2)^2-aa2(i)^2)^0.5));
   %滑滚比S22=deltaU2/U2，滚动方向
    S22(i)=abs(deltaU2(i)/U2(i));
    
    U2r(i)=U2(i)+0.5*U2(i)*S22(i);
    U2b(i)=U2(i)-0.5*U2(i)*S22(i);
    
    S2(i)=(S21(i)+S22(i))/2;
    P02(i)=Q2(i)/Rx2^2/E2;
    U02(i)=niandu2(i)*U2(i)/(E2*Rx2);
    sita2(i)=sita0*(K/niandu2(i)/(U2(i)^2));
%     A2(i)=-1.3527*10^(-20)*P02(i)^1.0849*U02(i)^(-0.28538)*sita2(i)^(-1.8235);%博士论文第四章
%     B2(i)=1.0217*10^(-25)*P02(i)^2.8202*U02(i)^(-0.53722)*sita2(i)^(-2.71);
%     C2(i)=1.6836*10^(-23)*P02(i)^1.7563*U02(i)^(-0.25749)*sita2(i)^(-2.6437);
%     D2(i)=8.0232*10^(-23)*P02(i)^1.1072*U02(i)^(-0.3333)*sita2(i)^(-2.0078);    
    
    A2(i)=-6.96*10^(2)*P02(i)^(-0.163)*U02(i)^(0.577)*sita2(i)^(0.308);%王艳霜4109发表论文
    B2(i)=5.19*10^(6)*P02(i)^(-0.373)*U02(i)^(1.000)*sita2(i)^(0.687);
    C2(i)=6.826*10^(5)*P02(i)^(-0.170)*U02(i)^(0.578)*sita2(i)^(0.373);
    D2(i)=6.96*10^(2)*P02(i)^(-0.163)*U02(i)^(0.577)*sita2(i)^(0.308); 
    
    %A2(i)=-0.0196;
    %B2(i)=0.39873;
    %C2(i)=39.002;
    %D2(i)=0.01984;
    
    miu2(i)=0.0127*(50/(50-S22(i)))*log(0.584*Q2(i)/niandu0/deltaU2(i)/(U2(i))^2);  % 滚动体/内圈的拖动系数!   只是用来计算托动力，不影响其他%CHENZHEGAI
    T2(i)=miu2(i)*Q2(i);     
else
    T1(i)=fcbr*Q1(i);
    T2(i)=fcbr*Q2(i);
    Re=2*f1*Dw/(2*f1+1);   
    deltaU1(i)=abs(-Wo(i)*Dm/2+(abs(Wx(i))*cos(a1(i))+ abs(Wz(i))*sin(a1(i))-Wo(i)*cos(a1(i)))*(Re-(Re^2-aa1(i)^2)^0.5+((Dw/2)^2-aa1(i)^2)^0.5));
    U1(i)=abs(-Wo(i)*Dm/4-0.5*(abs(Wx(i))*cos(a1(i))+ abs(Wz(i))*sin(a1(i))+Wo(i)*cos(a1(i)))*(Re-(Re^2-aa1(i)^2)^0.5+((Dw/2)^2-aa1(i)^2)^0.5));
    Ri=2*f2*Dw/(2*f2+1);
    deltaU2(i)=abs(+(W2-Wo(i))*Dm/2-( abs(Wx(i))*cos(a2(i))+ abs (Wz(i))*sin(a2(i))+(W2-Wo(i))*cos(a2(i)))*(Ri-(Ri^2-aa2(i)^2)^0.5+((Dw/2)^2-aa2(i)^2)^0.5));
    U2(i)=abs(+(W2-Wo(i))*Dm/4+0.5*( abs(Wx(i))*cos(a2(i))+abs(Wz(i))*sin(a2(i))-(W2-Wo(i))*cos(a2(i)))*(Ri-(Ri^2-aa2(i)^2)^0.5+((Dw/2)^2-aa2(i)^2)^0.5));    
end
% 最小油膜厚度计算!!!
L1(i)=niandu0*(U1(i))^2*beita0/K;       
Ct1(i)=1/(1+0.241*(1+14.8*S1(i)^0.83)*(L1(i))^0.64);  %Chenzhe Gai 热修正因子
L2(i)=niandu0*(U2(i))^2*beita0/K;
Ct2(i)=1/(1+0.241*(1+14.8*S2(i)^0.83)*(L2(i))^0.64);micro_config = load_micro_interface_config();
thermal_factor1 = 1; thermal_factor2 = 1; texture_factor = 1;
if micro_config.thermal.enabled
    thermal_factor1 = micro_config.thermal.Ct1;
    thermal_factor2 = micro_config.thermal.Ct2;
end
if micro_config.texture.enabled
    texture_factor = micro_config.texture.Cr;
end
oilh1(i)=thermal_factor1*texture_factor*Rx1*3.63*(niandu0*U1(i)/(E1*Rx1))^(0.68)*(nianya0*E1)^(0.49)*(Q1(i)/Rx1^2/E1)^(-0.073)*(1-exp(-0.68*K1));
oilh2(i)=thermal_factor2*texture_factor*Rx2*3.63*(niandu0*U2(i)/(E2*Rx2))^(0.68)*(nianya0*E2)^(0.49)*(Q2(i)/Rx2^2/E2)^(-0.073)*(1-exp(-0.68*K2));    

miuI(i)=0.0127*(50/(50-S12(i)))*log(0.584*Q1(i)/niandu0/deltaU1(i)/(U1(i))^2); %CHENZHEGAI   Bennit公式 
miuO(i)=0.0127*(50/(50-S22(i)))*log(0.584*Q2(i)/niandu0/deltaU2(i)/(U2(i))^2);
fs1(i)=miuI(i)*Q1(i);  fs2(i)=miuO(i)*Q2(i);

Rbp=2/(2/Dw-1/Rp);
%求球/保持架之间的法向力和摩擦力！
Kbp=11/(Rp-Dw/2);
deltacx(i)=pi*Dm/n*(Wo(i)/Wc-1);     % 滚动体偏离原来位置的距离！
Fsrx(i)=Kbp*deltacx(i);
fstx(i)=0.05*Fsrx(i);

if niandu0==0
    Fsrz(i)=0;
    fstz(i)=0;
else
    deltacz(i)=deltac0-sqrt((x(i)-0)^2);  deltacz(i)=abs(deltacz(i));
    Fsrz(i)=2.245*pi*niandu0*abs(Wz(i))*Dw/2*Rbp^0.5*Dw/2*deltacz(i)^(-0.5)*(x(i)-0)/abs(x(i)-0);
    fstz(i)=7.3*niandu0*abs(Wz(i))*Dw/2*sqrt(Rbp*Dw/2)*log(4*abs(Fsrz(i))/(4.49*pi*niandu0*abs(Wz(i))*Dw/2*Rbp))*(x(i)-0)/abs(x(i)-0);
end

%阻力损失计算！
if lubricationtype==0
    Ree(i)=oilden*Dw/niandu01(i)*0.5*Dm*Wo(i) ;                        %oilden为油的密度！Ree是估计的雷诺数！
else
    if niandu0==0
        Ree(i)=0;
    else
        Ree(i)=oilden*Dw/niandu0*0.5*Dm*Wo(i) ;                        %oilden为油的密度！Ree是估计的雷诺数！
    end
end
if Ree(i)==0
    Cd(i)=0;
elseif Ree(i)<=0.1
    Cd(i)=275;
elseif Ree(i)<=1
    Cd(i)=interp1([-1 0],[275 30],log10(Ree(i)));
elseif Ree(i)<=10
    Cd(i)=interp1([0 1],[30 4.2],log10(Ree(i)));
elseif Ree(i)<=100
    Cd(i)=interp1([1 2],[4.2 1.2],log10(Ree(i)));
elseif Ree(i)<=1000
    Cd(i)=interp1([2 3],[1.2 0.48],log10(Ree(i)));
elseif Ree(i)<=10000
    Cd(i)=interp1([3 4],[0.48 0.4],log10(Ree(i)));
elseif Ree(i)<=100000
    Cd(i)=interp1([4 5],[0.4 0.45],log10(Ree(i)));
elseif Ree(i)<=200000
    Cd(i)=interp1([5 5.301],[0.45 0.4],log10(Ree(i)));
elseif Ree(i)<=300000
    Cd(i)=interp1([5.301 5.47710],[0.4 0.1],log10(Ree(i)));
elseif Ree(i)<=400000
    Cd(i)=interp1([5.47710 5.60210],[0.1 0.09],log10(Ree(i)));
elseif Ree(i)<=500000
    Cd(i)=interp1([5.60210 5.69900],[0.09 0.09],log10(Ree(i)));
elseif Ree(i)<=1000000
    Cd(i)=interp1([5.69900 6],[0.09 0.09],log10(Ree(i)));
else
    Cd(i)=0.09;
end
Fzu(i)=0.015*oilden*pi*Cd(i)*Dw^2*(Wo(i)*Dm)^1.95/10;      %求出滚动体受到的阻力！

 %速度求解方程包括由转速引起的求惯性力和惯性力矩项
if loadj==1
Fy(i)=m*(r(i)+Dm/2)*(Wononload-Wo(i))/(2*pi/n)*Wo(i);   %求惯性力！
z2(i)=Fsrx(i)+Fzu(i)+T1(i)-T2(i)+Fy(i)-fstz(i);%球平衡方程（b）
z4(i)=-T1(i)*Dw/2*cos(a1(i))-T2(i)*Dw/2*cos(a2(i))+fstx(i)*Dw/2+J*(Wxnonload-Wx(i))/(2*pi/n)*Wo(i); %球平衡方程（d，e）   
z5(i)=-(fs1(i)+fs2(i))*Dw/2*Wy(i)/abs(Wy(i))-(J*Wo(i)*Wz(i));%球平衡方程（f）   
else
if i==loadi(1)
Fy(i)=m*(r(i)+Dm/2)*(Wo(i+1)-Wo(i))/(2*pi/n)*Wo(i);   %求钢球的惯性力！
z2(i)=Fsrx(i)+Fzu(i)+T1(i)-T2(i)+Fy(i)-fstz(i);%球平衡方程（b）
z4(i)=-T1(i)*Dw/2*cos(a1(i))-T2(i)*Dw/2*cos(a2(i))+fstx(i)*Dw/2+J*(Wx(i+1)-Wx(i))/(2*pi/n)*Wo(i);  %球平衡方程（d，e）  
z5(i)=-(fs1(i)+fs2(i))*Dw/2*Wy(i)/abs(Wy(i))-(J*Wo(i)*Wz(i))+J*(Wy(i+1)-Wy(i))/(2*pi/n)*Wo(i); %球平衡方程（f） 
else
    if Q1(i)>Q1(i-1)
      Fy(i)=m*(r(i)+Dm/2)*(Wo(i-1)-Wo(i))/( 2*pi/n)*Wo(i); %求钢球的惯性力！
      z2(i)=Fsrx(i)+Fzu(i)+T1(i)-T2(i)+Fy(i)-fstz(i);
      z4(i)=-T1(i)*Dw/2*cos(a1(i))-T2(i)*Dw/2*cos(a2(i))+fstx(i)*Dw/2-J*(Wx(i-1)-Wx(i))/( 2*pi/n)*Wo(i);
      z5(i)=-(fs1(i)+fs2(i))*Dw/2*Wy(i)/abs(Wy(i))-(J*Wo(i)*Wz(i))-J*(Wy(i-1)-Wy(i))/( 2*pi/n)*Wo(i);
    else
      Fy(i)=m*(r(i)+Dm/2)*(Wo(i)-Wo(i-1))/( 2*pi/n)*Wo(i); %求钢球的惯性力！
      z2(i)=Fsrx(i)+Fzu(i)+T1(i)-T2(i)+Fy(i)-fstz(i);
      z4(i)=-T1(i)*Dw/2*cos(a1(i))-T2(i)*Dw/2*cos(a2(i))+fstx(i)*Dw/2-J*(Wx(i)-Wx(i-1))/( 2*pi/n)*Wo(i);
      z5(i)=-(fs1(i)+fs2(i))*Dw/2*Wy(i)/abs(Wy(i))-(J*Wo(i)*Wz(i))-J*(Wy(i)-Wy(i-1))/( 2*pi/n)*Wo(i);
    end
end
end
end
  
if loadj==n                           %  此时没有非承载区！！
    z2nonload=0;
    z4nonload=0;
    pvzhi1nonload=0;
else                                   %  存在没有非承载区！！
% 计算非承载区滚子的转速（认为非承载区各滚动体运动状态相同！）
%  非承载区载荷计算
a1nonload=0;    % 与外圈接触角为0
Q1nonload=m*Wononload^2*Dm/2;  fs1nonload=0.02*Q1nonload;  save Q1sortnonload;
aa1nonload=(6*K1^2*ee21*R21*Q1nonload/(E1*pi))^0.3333; b1nonload=aa1nonload/K1;    %求得接触椭圆长短轴！
Ph1nonload=1.5*Q1nonload/(pi*aa1nonload*b1nonload);   save Ph1sortnonload
if lubricationtype==0
    %黏度计算
    niandu00=0.132*sinh(925.19/(-148.85+sita0+273))*oilden/1e6;      %120度时的粘度，wys第40页
    niandu100=niandu00*exp(Ph1nonload*nianya0);      % 考虑粘压系数后的 黏度！
    %求拖动系数
    Re=2*f1*Dw/(2*f1+1);   
    deltaU1nonload=abs(-Wononload*Dm/2+(abs(Wxnonload)*cos(a1nonload)+ abs(Wznonload)*sin(a1nonload)-Wononload*cos(a1nonload))*((Re^2-0*aa1nonload^2)^0.5-(Re^2-aa1nonload^2)^0.5+((Dw/2)^2-aa1nonload^2)^0.5));
    U1nonload=abs(-Wononload*Dm/4-0.5*(abs(Wxnonload)*cos(a1nonload)+ abs(Wznonload)*sin(a1nonload)+Wononload*cos(a1nonload))*((Re^2-0*aa1nonload^2)^0.5-(Re^2-aa1nonload^2)^0.5+((Dw/2)^2-aa1nonload^2)^0.5));
    S11nonload=deltaU1nonload/U1nonload;
    deltaU1nonload=abs(-Wononload*Dm/2+(abs(Wxnonload)*cos(a1nonload)+ abs(Wznonload)*sin(a1nonload)-Wononload*cos(a1nonload))*((Re^2-aa1nonload^2)^0.5-(Re^2-aa1nonload^2)^0.5+((Dw/2)^2-aa1nonload^2)^0.5));
    U1nonload=abs(-Wononload*Dm/4-0.5*(abs(Wxnonload)*cos(a1nonload)+ abs(Wznonload)*sin(a1nonload)+Wononload*cos(a1nonload))*((Re^2-aa1nonload^2)^0.5-(Re^2-aa1nonload^2)^0.5+((Dw/2)^2-aa1nonload^2)^0.5));
    S12nonload=abs(deltaU1nonload/U1nonload);
    S1nonload=(S11nonload+S12nonload)/2;
    P01nonload=Ph1nonload/E1;
    U01nonload=niandu100*U1nonload/(E1*Rx1);
    sita1nonload=sita0*sqrt(K*beita0*niandu100)/(E1*Rx1);
    A1nonload=-1.3527*10^(-22)*P01nonload^1.0849*U01nonload^(-0.28538)*sita1nonload^(-1.8235);
    B1nonload=1.0217*10^(-25)*P01nonload^2.8202*U01nonload^(-0.53722)*sita1nonload^(-2.71);
    C1nonload=1.6836*10^(-23)*P01nonload^1.7563*U01nonload^(-0.25749)*sita1nonload^(-2.6437);
    D1nonload=8.0232*10^(-23)*P01nonload^1.1072*U01nonload^(-0.3333)*sita1nonload^(-2.0078);
    miu1nonload=(A1nonload+B1nonload*S1nonload)*exp(-C1nonload*S1nonload)+D1nonload;  % 球/外圈的拖动系数!
    T1nonload=miu1nonload*Q1nonload;                    % 球/外圈的拖动力！
else
    Re=2*f1*Dw/(2*f1+1);   
    deltaU1nonload=abs(-Wononload*Dm/2+(abs(Wxnonload)*cos(a1nonload)+ abs(Wznonload)*sin(a1nonload)-Wononload*cos(a1nonload))*((Re^2-0*aa1nonload^2)^0.5-(Re^2-aa1nonload^2)^0.5+((Dw/2)^2-aa1nonload^2)^0.5));
    T1nonload=fcbr*Q1nonload;
end   
    
%求球/保持架之间的法向力和摩擦力！
Kbp=11/(Rp-Dw/2);
deltacxnonload=pi*Dm/n*(1-Wononload/Wc);   % 球偏离原来位置的距离！
Fsrxnonload=Kbp*deltacxnonload;    fstxnonload=0.05*Fsrxnonload;
Fsrznonload=0;   fstznonload=0;
%阻力损失计算！
if lubricationtype==0
    Reenonload=oilden*Dw/niandu00*0.5*Dm*Wononload ;                 %oilden为油的密度！Ree是估计的雷诺数！
else
    if niandu0==0
        Reenonload=0;
    else
        Reenonload=oilden*Dw/niandu0*0.5*Dm*Wononload ;                 %oilden为油的密度！Ree是估计的雷诺数！
    end
end
%插值求阻力系数！
if Reenonload==0
    Cdnonload=0;
elseif Reenonload<=0.1
    Cdnonload=275;
elseif Reenonload<=1
    Cdnonload=interp1([-1 0],[275 30],log10(Reenonload));
elseif Reenonload<=10
    Cdnonload=interp1([0 1],[30 4.2],log10(Reenonload));
elseif Reenonload<=100
    Cdnonload=interp1([1 2],[4.2 1.2],log10(Reenonload));
elseif Reenonload<=1000
    Cdnonload=interp1([2 3],[1.2 0.48],log10(Reenonload));
elseif Reenonload<=10000
    Cdnonload=interp1([3 4],[0.48 0.4],log10(Reenonload));
elseif Reenonload<=100000
    Cdnonload=interp1([4 5],[0.4 0.45],log10(Reenonload));
elseif Reenonload<=200000
    Cdnonload=interp1([5 5.301],[0.45 0.4],log10(Reenonload));
elseif Reenonload<=300000
    Cdnonload=interp1([5.301 5.47710],[0.4 0.1],log10(Reenonload));
elseif Reenonload<=400000
    Cdnonload=interp1([5.47710 5.60210],[0.1 0.09],log10(Reenonload));
elseif Reenonload<=500000
    Cdnonload=interp1([5.60210 5.69900],[0.09 0.09],log10(Reenonload));
elseif Reenonload<=1000000
    Cdnonload=interp1([5.69900 6],[0.09 0.09],log10(Reenonload));
else
    Cdnonload=0.09;
end
Fzunonload=0.015*oilden*pi*Cdnonload*Dw^2*(Wononload*Dm)^1.95/10;      % 0.015为考虑油气混合物后添加的系数  求出滚动体受到的阻力！
z2nonload=-Fsrxnonload+Fzunonload+T1nonload;
z4nonload=-T1nonload*Dw/2*cos(a1nonload)+fstxnonload*Dw/2;
pvzhi1nonload=Ph1nonload*deltaU1nonload; 
end

zzzz34=[z2';z4';z5';z2nonload*(n-loadj)*10;z4nonload*(n-loadj)*10;];
save zzzz34;
nnnn=length(zzzz34);

result1134=0;
for i=1:nnnn
    result1134=result1134+zzzz34(i)^2;
end
save result1134
save T1;save T2;save Fsrx;save fstx;save S1; save S2

for i=1:loadj
pvzhi1(i)=Ph1(i)*deltaU1(i);  pvzhi2(i)=Ph2(i)*deltaU2(i);
end
save pvzhi1; save pvzhi2; save pvzhi1nonload;save deltaU1; save deltaU2;
save oilh1;save oilh2;                 %最小油膜厚度！！






