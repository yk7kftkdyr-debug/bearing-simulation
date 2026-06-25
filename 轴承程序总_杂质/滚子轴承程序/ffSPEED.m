function ffSPEED(datafromvb,ttt)     % 求Wo Wx Wc=mean(Wo)

rollerspeed=ttt;
load Ph1ii;load Ph2ii;load Q1ii;load Q2ii;
load loadii;

for i=1:loadj
    if Q1ii(i)==max(Q1ii)
        markQ1=i;          % 找到受载最大的滚子
    end
end

% 数据输入部分！！
n=datafromvb(1);Dw=datafromvb(2);Dm=datafromvb(3);lenroller=datafromvb(4);lenrollerline=datafromvb(5);arcroller=datafromvb(6);
deltar0=datafromvb(7); Rp=datafromvb(8); %兜孔半径
tuoyuandu=datafromvb(9); %外圈的椭圆度
ballden=datafromvb(59);e1=datafromvb(11);e2=datafromvb(12);e3=datafromvb(13);o1=datafromvb(14);o2=datafromvb(15);o3=datafromvb(16);
W1=datafromvb(17)*pi/30;W2=datafromvb(18)*pi/30; Fyy=datafromvb(19); Fzz=datafromvb(20)+0.1; Myy=datafromvb(21)+0.1; Mzz=datafromvb(22)+0.1;
cucao1=datafromvb(23);cucao2=datafromvb(24);cucao3=datafromvb(25);
%润滑油基本参数
oilden=datafromvb(26); sita0=datafromvb(27);K=datafromvb(28); %常温下的导热系数
niandu0=datafromvb(29); nianya0=datafromvb(30); beita0=datafromvb(31);   %粘温系数（近似认为不变，wys论文中）   所有参数均为常温下的参数！！！
yindao=datafromvb(32); yindaojianxi=datafromvb(33);   % 引导方式：引导间隙
lubricationtype=datafromvb(70);
fcbr=datafromvb(71); fcbc=datafromvb(72); fccr=datafromvb(73);
Dr1=Dm+Dw+deltar0;Dr2=Dm-Dw-deltar0;  % 滚道直径
gama=Dw/Dm;

Rr1=Dr1/2*Dw/2/(Dr1/2-Dw/2);  % 滚子与外圈的当量曲率半径！！！
Rr2=Dr2/2*Dw/2/(Dr2/2+Dw/2);  % 滚子与内圈的当量曲率半径！！！
E1=2/(((1-o1^2)/e1)+((1-o3^2)/e3));             %e1，e3为外圈和球的弹性模量，o1，o3为外圈，球的泊松比
E2=2/(((1-o2^2)/e2)+((1-o3^2)/e3));             %e1，e3为内圈和球的弹性模量，o1，o3为内圈，球的泊松比

m=ballden*Dw*Dw*lenroller*pi/4;  %滚子的质量计算
J=m*Dw*Dw/8;             %转动惯量计算
for i=1:loadj
    Wo(i)=rollerspeed(i);      Wx(i)=rollerspeed(loadj+i);
end
Wononload=rollerspeed(2*loadj+1); Wxnonload=rollerspeed(2*loadj+2);
Wc=(sum(Wo)+Wononload*(n-loadj))/n;          % 保持架转速取所有滚子转速的平均值！！


for i=1:loadj       %共n个球循环
    sita(i)=2*pi*(loadii(i))/n;
end
for i=1:loadj
    if lubricationtype==0
        % 黏度计算！！
        oilden1(i)=oilden*(1+(0.6e-9*Ph1ii(i))/(1+1.7e-9*Ph1ii(i)) );   % 滑油密度计算    公式在(清华)弹流12页
        niandu01(i)=0.132*sinh(925.19/(-148.85+sita0+273))*oilden1(i)/1e6;          %sita0 度时的粘度，wys第40页
        niandu1(i)=niandu01(i)*exp(Ph1ii(i)*nianya0);      % 考虑粘压系数后的 黏度！

        oilden2(i)=oilden*(1+(0.6e-9*Ph2ii(i))/(1+1.7e-9*Ph2ii(i)) );   % 滑油密度计算    公式在清华弹流12页
        niandu02(i)=0.132*sinh(925.19/(-148.85+sita0+273))*oilden2(i)/1e6; % sita0度时的粘度，wys第40页
        niandu2(i)=niandu02(i)*exp(Ph2ii(i)*nianya0);      % 考虑粘压系数后的 黏度！

        %求拖动系数
        deltaU1(i)=abs(Dm/2*((1+gama)*(W1+Wo(i))-gama*Wx(i)));
        U1(i)=0.5*0.5*Dm*((1+gama)*(W1+Wo(i))+gama*Wx(i));
        S1(i)=deltaU1(i)/U1(i);                 % 滑滚比
        P01(i)=Ph1ii(i)/E1;
        U01(i)=niandu1(i)*U1(i)/(E1*Rr1);
        sita1(i)=sita0*sqrt(K*beita0*niandu1(i))/(E1*Rr1);
        A1(i)=-1.3527*10^(-22)*P01(i)^1.0849*U01(i)^(-0.28538)*sita1(i)^(-1.8235);
        B1(i)=1.0217*10^(-25)*P01(i)^2.8202*U01(i)^(-0.53722)*sita1(i)^(-2.71);
        C1(i)=1.6836*10^(-23)*P01(i)^1.7563*U01(i)^(-0.25749)*sita1(i)^(-2.6437);
        D1(i)=8.0232*10^(-23)*P01(i)^1.1072*U01(i)^(-0.3333)*sita1(i)^(-2.0078);
        miu1(i)=(A1(i)+B1(i)*S1(i))*exp(-C1(i)*S1(i))+D1(i);  % 滚子/外圈的拖动系数!
        T1(i)=miu1(i)*Q1ii(i);                    % 滚子/外圈的拖动力！

        deltaU2(i)=abs(Dm/2*((1-gama)*(W2-Wo(i))-gama*Wx(i)));
        U2(i)=0.5*0.5*Dm*((1-gama)*(W2-Wo(i))+gama*Wx(i));
        S2(i)=deltaU2(i)/U2(i);                 % 滑滚比
        P02(i)=Ph2ii(i)/E2;
        U02(i)=niandu2(i)*U2(i)/(E2*Rr2);
        sita2(i)=sita0*sqrt(K*beita0*niandu2(i))/(E2*Rr2);
        A2(i)=-1.3527*10^(-22)*P02(i)^1.0849*U02(i)^(-0.28538)*sita2(i)^(-1.8235);
        B2(i)=1.0217*10^(-25)*P02(i)^2.8202*U02(i)^(-0.53722)*sita2(i)^(-2.71);
        C2(i)=1.6836*10^(-23)*P02(i)^1.7563*U02(i)^(-0.25749)*sita2(i)^(-2.6437);
        D2(i)=8.0232*10^(-23)*P02(i)^1.1072*U02(i)^(-0.3333)*sita2(i)^(-2.0078);
        miu2(i)=(A2(i)+B2(i)*S2(i))*exp(-C2(i)*S2(i))+D2(i);  % 滚子/内圈的拖动系数!
        T2(i)=miu2(i)*Q2ii(i);
    else
        T1(i)=fcbr*Q1ii(i);
        T2(i)=fcbr*Q2ii(i);
        U1(i)=0.5*0.5*Dm*((1+gama)*(W1+Wo(i))+gama*Wx(i));
        U2(i)=0.5*0.5*Dm*((1-gama)*(W2-Wo(i))+gama*Wx(i));
        deltaU1(i)=abs(Dm/2*((1+gama)*(W1+Wo(i))-gama*Wx(i)));
        deltaU2(i)=abs(Dm/2*((1-gama)*(W2-Wo(i))-gama*Wx(i)));
    end
end

   % 最小油膜厚度计算!!!
   %yang解释：等温、考虑热效应、考虑表面纹理参数，使用对应油膜厚度计算公式；
F_radial = sqrt(Fyy^2 + Fzz^2);
thetaF = atan2(Fyy, Fzz);

Q20 = zeros(1, loadj);
Q10 = zeros(1, loadj);
Fc0 = zeros(1, loadj);

for i = 1:loadj
    sita(i) = 2*pi*loadi(i)/n;
    thetaRel = atan2(sin(sita(i)-thetaF), cos(sita(i)-thetaF));

    Qmax = 4*F_radial/(pi*n);

    if cos(thetaRel) > 0
        Q20(i) = Qmax*cos(thetaRel);
    else
        Q20(i) = 0;
    end

    Fc0(i) = m*Wo(i)^2*Dm/2;
    Q10(i) = Q20(i) + Fc0(i);

    if Q20(i) <= 0
        Q20(i) = 1;
    end

    if Q10(i) <= 0
        Q10(i) = 1;
    end
end
%载荷分布！！( 假设承载区载荷对称分布！！！) 
thickness=lenroller/150;
UU1=0.5*0.5*Dm*((1+gama)*(W1+(W2*(1-gama)-W1*(1+gama))/2)+gama*(W2-(W2*(1-gama)-W1*(1+gama))/2)*(1-gama)/gama);
UU2=0.5*0.5*Dm*((1-gama)*(W2-(W2*(1-gama)-W1*(1+gama))/2)+gama*(W2-(W2*(1-gama)-W1*(1+gama))/2)*(1-gama)/gama);

widthcontact11=zeros(150,loadj);  widthcontact22=zeros(150,loadj);      % 定义接触宽度矩阵！！！
Ph11=zeros(150,loadj);  Ph22=zeros(150,loadj);      % 定义接触应力矩阵！！！

for i=1:loadj              %包括了0号滚子，
     sita(i)=2*pi*(loadi(i))/n;   
     Q1e(i)=0;Q2e(i)=0; FP1e(i)=0;FP2e(i)=0; eQ2fenzi(i)=0; eQ2fenmu(i)=0;  
     numcontact11(i)=0;  numcontact22(i)=0;     %   统计承载的片数！！
     
   % 最小油膜厚度计算!!!
   %yang解释：等温、考虑热效应、考虑表面纹理参数，使用对应油膜厚度计算公式；考虑杂质影响，在接触变形delta1k，delta2k上加上特征位移ud  
   
%***************************等温条件下膜厚计算公式******************************************      
 U1film=max(abs(U1(i)),eps); U2film=max(abs(U2(i)),eps);
 Q1film=max(Q1ii(i),1); Q2film=max(Q2ii(i),1);micro_config = load_micro_interface_config();
thermal_factor1 = 1; thermal_factor2 = 1; texture_factor = 1; debris_shift = 0;
if micro_config.thermal.enabled
    thermal_factor1 = micro_config.thermal.Ct1;
    thermal_factor2 = micro_config.thermal.Ct2;
end
if micro_config.texture.enabled
    texture_factor = micro_config.texture.Cr;
end
if micro_config.debris.enabled
    debris_shift = micro_config.debris.ud + micro_config.debris.debris_displacement;
end
oilh1(i)=thermal_factor1*texture_factor*2.65*Rr1*(nianya0*E1)^0.54*(niandu0*U1film/(E1*Rr1))^0.7*(Q1film/(lenroller*E1*Rr1))^(-0.13);oilh2(i)=thermal_factor2*texture_factor*2.65*Rr2*(nianya0*E2)^0.54*(niandu0*U2film/(E2*Rr2))^0.7*(Q2film/(lenroller*E2*Rr2))^(-0.13);
 
%***********************************考虑热效应膜厚计算公式******************************************************   
% thermal correction option is controlled by micro_config.thermal.
% thermal correction option is controlled by micro_config.thermal.
% optional thermal film correction uses micro_config.thermal.
% optional thermal film correction uses micro_config.thermal.
%*****************************************************************************************   

%************************************考虑表面纹理参数膜厚计算公式*****************************************************   
% texture option is controlled by micro_config.texture.
% texture correction factor is controlled by micro_config.texture.
% optional texture film correction uses micro_config.texture.
% optional texture film correction uses micro_config.texture.
% %*******************************************************************************************************   
%*******************************考虑杂质影响的特征位移**********************************************************   
% debris displacement option is controlled by micro_config.debris.
%oilh1(i)=2.65*nianya0^0.54*(niandu0*UU1)^0.7*Rr1^0.43/((Q10(i)/lenroller/E1/Rr1)^0.13*E1^0.03);
%oilh2(i)=2.65*nianya0^0.54*(niandu0*UU2)^0.7*Rr2^0.43/((Q20(i)/lenroller/E2/Rr2)^0.13*E2^0.03);
%*****************************************************************************************   
    %求滚子/保持架之间的法向力和摩擦力！
    Kbp=11/(Rp-Dw/2)/3;
    deltacx(i)=pi*Dm/n*(Wo(i)/Wc-1);     % 滚子偏离原来位置的距离！
    Fm(i)=Kbp*(deltacx(i));
    fm(i)=0.05*Fm(i);

    %阻力损失计算！
    if lubricationtype==0
        Ree(i)=oilden*Dw/niandu01(i)*0.5*Dm*Wo(i) ;    %oilden为油的密度！Ree是估计的雷诺数！
    else
        Ree(i)=oilden*Dw/niandu0*0.5*Dm*Wo(i) ;
    end
    %插值求阻力系数！
    if Ree(i)<=0.1
        Cd(i)=60;
    elseif Ree(i)<=1
        Cd(i)=interp1([-1 0],[60 10],log10(Ree(i)));
    elseif Ree(i)<=10
        Cd(i)=interp1([0 1],[10 3],log10(Ree(i)));
    elseif Ree(i)<=100
        Cd(i)=interp1([1 2],[3 1.8],log10(Ree(i)));
    elseif Ree(i)<=1000
        Cd(i)=interp1([2 3],[1.8 1],log10(Ree(i)));
    elseif Ree(i)<=10000
        Cd(i)=interp1([3 4],[1 1.2],log10(Ree(i)));
    elseif Ree(i)<=200000
        Cd(i)=interp1([4 5.301],[1.2 1.2],log10(Ree(i)));
    elseif Ree(i)<=300000
        Cd(i)=interp1([5.301 5.47710],[1.2 0.9],log10(Ree(i)));
    elseif Ree(i)<=400000
        Cd(i)=interp1([5.47710 5.60210],[0.9 0.65],log10(Ree(i)));
    elseif Ree(i)<=500000
        Cd(i)=interp1([5.60210 5.69900],[0.65 0.3],log10(Ree(i)));
    elseif Ree(i)<=1000000
        Cd(i)=interp1([5.69900 6],[0.3 0.3],log10(Ree(i)));
    else
        Cd(i)=0.3;
    end
    Fzu(i)=0.015*oilden*Cd(i)*(Wo(i)*Dm/2)^2*lenroller*(Dw-Dw/2);      % 0.015为考虑油气混合物后添加的系数  求出滚动体受到的阻力！

    %  动压力计算！！
    if lubricationtype==0
        FP1(i)=lenroller*E1*Rr1*1.43*(1+gama)*(niandu01(i)*U1(i)/(2*E1*Rr1))^0.71;
        FP2(i)=lenroller*E2*Rr2*1.43*(1+gama)*(niandu02(i)*U2(i)/(2*E2*Rr2))^0.71;
    else
        FP1(i)=lenroller*E1*Rr1*1.43*(1+gama)*(niandu0*U1(i)/(2*E1*Rr1))^0.71;
        FP2(i)=lenroller*E2*Rr2*1.43*(1+gama)*(niandu0*U2(i)/(2*E2*Rr2))^0.71;
    end


    if loadj==1
        Fy(i)=m*Dm/2*(Wononload-Wo(i))/(2*pi/n)*Wo(i);   %求惯性力！
        z1(i)=T1(i)+FP1(i)-FP2(i)-T2(i)+Fm(i)+Fy(i)+Fzu(i);
        z3(i)=Dw/2*(T2(i)+T1(i)-fm(i))-J*(Wxnonload-Wx(i))/(2*pi/n)*Wo(i);
    else
        if i==1
            Fy(i)=m*Dm/2*(Wo(i+1)-Wo(i))/(2*pi/n)*Wo(i);   %求惯性力！
            z1(i)=T1(i)+FP1(i)-FP2(i)-T2(i)+Fm(i)+Fy(i)+Fzu(i);
            z3(i)=Dw/2*(T2(i)+T1(i)-fm(i))-J*(Wx(i+1)-Wx(i))/(2*pi/n)*Wo(i);
        else
            if i>=markQ1
                Fy(i)=m*Dm/2*(Wo(i)-Wo(i-1))/(2*pi/n)*Wo(i);   %求惯性力！
                z1(i)=T1(i)+FP1(i)-FP2(i)-T2(i)+Fm(i)+Fy(i)+Fzu(i);
                z3(i)=Dw/2*(T2(i)+T1(i)-fm(i))-J*(Wx(i)-Wx(i-1))/(2*pi/n)*Wo(i);
            else
                Fy(i)=m*Dm/2*(Wo(i-1)-Wo(i))/(2*pi/n)*Wo(i);   %求惯性力！
                z1(i)=T1(i)+FP1(i)-FP2(i)-T2(i)+Fm(i)+Fy(i)+Fzu(i);
                z3(i)=Dw/2*(T2(i)+T1(i)-fm(i))-J*(Wx(i-1)-Wx(i))/(2*pi/n)*Wo(i);
            end
        end
    end
end  % 最后一个end, 对应于最外层的for循环！！


% 计算非承载区滚子的转速（认为非承载区各滚子运动状态相同！）
%  非承载区载荷计算
Q1nonload=0.5*m*Dm*Wononload^2;  save Q1nonload;
widthcontactnonload=sqrt(8*Q1nonload*Rr1/(pi*lenroller*E1));  save widthcontactnonload;
Ph1nonload=sqrt(Q1nonload*E1/(2*pi*lenroller*Rr1));  save Ph1nonload
if lubricationtype==0
    %黏度计算
    oilden00=oilden*(1+(0.6e-9*Ph1nonload)/(1+1.7e-9*Ph1nonload) );   % 滑油密度计算    公式在清华弹流12页
    niandu00=0.132*sinh(925.19/(-148.85+sita0+273))*oilden00/1e6;     %sita0 度时的粘度，wys第40页
    niandu100=niandu00*exp(Ph1nonload*nianya0);      % 考虑粘压系数后的 黏度！

    % 求解拖动力
    deltaU1nonload=abs(Dm/2*((1+gama)*(W1+Wononload)-gama*Wxnonload));
    U1nonload=0.5*0.5*Dm*((1+gama)*(W1+Wononload)+gama*Wxnonload);
    S1nonload=deltaU1nonload/U1nonload ;               % 滑滚比
    P01nonload=Ph1nonload/E1;
    U01nonload=niandu100*U1nonload/(E1*Rr1);
    sita1nonload=sita0*sqrt(K*beita0*niandu100)/(E1*Rr1);
    A1nonload=-1.3527*10^(-22)*P01nonload^1.0849*U01nonload^(-0.28538)*sita1nonload^(-1.8235);
    B1nonload=1.0217*10^(-25)*P01nonload^2.8202*U01nonload^(-0.53722)*sita1nonload^(-2.7);
    C1nonload=1.6836*10^(-23)*P01nonload^1.7563*U01nonload^(-0.25749)*sita1nonload^(-2.6437);
    D1nonload=8.0232*10^(-23)*P01nonload^1.1072*U01nonload^(-0.3333)*sita1nonload^(-2.0078);
    miu1nonload=(A1nonload+B1nonload*S1nonload)*exp(-C1nonload*S1nonload)+D1nonload;       % 滚子/外圈的拖动系数!
    T1nonload=miu1nonload*Q1nonload  ;                % 滚子/外圈的拖动力！
else
    T1nonload=fcbr*Q1nonload  ;                % 滚子/外圈的拖动力！
    deltaU1nonload=abs(Dm/2*((1+gama)*(W1+Wononload)-gama*Wxnonload));
    U1nonload=0.5*0.5*Dm*((1+gama)*(W1+Wononload)+gama*Wxnonload);
end

%求球/保持架之间的法向力和摩擦力！
Kbp=11/(Rp-Dw/2)/3;
deltacxnonload=pi*Dm/n*(1-Wononload/Wc);     % 滚子偏离原来位置的距离！
Fmnonload=Kbp*deltacxnonload;
fmnonload=0.05*Fmnonload;

if lubricationtype==0
    FP1nonload=lenroller*E1*Rr1*1.43*(1+gama)*(niandu00*U1nonload/(2*E1*Rr1))^0.71;
else    
    FP1nonload=lenroller*E1*Rr1*1.43*(1+gama)*(niandu0*U1nonload/(2*E1*Rr1))^0.71;
end


%阻力损失计算！
if lubricationtype==0
    Reenonload=oilden*Dw/niandu00*0.5*Dm*Wononload ;     %oilden为油的密度！Ree是估计的雷诺数！
else
    Reenonload=oilden*Dw/niandu0*0.5*Dm*Wononload ;
end
%插值求阻力系数！
if Reenonload<=0.1
    Cdnonload=60;
elseif Reenonload<=1
    Cdnonload=interp1([-1 0],[60 10],log10(Reenonload));
elseif Reenonload<=10
    Cdnonload=interp1([0 1],[10 3],log10(Reenonload));
elseif Reenonload<=100
    Cdnonload=interp1([1 2],[3 1.8],log10(Reenonload));
elseif Reenonload<=1000
    Cdnonload=interp1([2 3],[1.8 1],log10(Reenonload));
elseif Reenonload<=10000
    Cdnonload=interp1([3 4],[1 1.2],log10(Reenonload));
elseif Reenonload<=200000
    Cdnonload=interp1([4 5.301],[1.2 1.2],log10(Reenonload));
elseif Reenonload<=300000
    Cdnonload=interp1([5.301 5.47710],[1.2 0.9],log10(Reenonload));
elseif Reenonload<=400000
    Cdnonload=interp1([5.47710 5.60210],[0.9 0.65],log10(Reenonload));
elseif Reenonload<=500000
    Cdnonload=interp1([5.60210 5.69900],[0.65 0.3],log10(Reenonload));
elseif Reenonload<=1000000
    Cdnonload=interp1([5.69900 6],[0.3 0.3],log10(Reenonload));
else
    Cdnonload=0.3;
end
Fzunonload=0.015*oilden*Cdnonload*(Wononload*Dm/2)^2*lenroller*(Dw-Dw/2) ;    % 0.015为考虑油气混合物后添加的系数  求出滚动体受到的阻力！
z1nonload=T1nonload+FP1nonload+Fzunonload -Fmnonload;
z3nonload=Dw/2*(T1nonload-fmnonload)-J*(Wxnonload-Wx(loadj))/(2*pi/n)*Wononload;

zzSPEED=[z1';z3';z1nonload*(n-loadj)*10;z3nonload*(n-loadj)*10];
save zzSPEED;
nnnn=length(zzSPEED);  result222=0;
for i=1:nnnn
    result222=result222+zzSPEED(i)^2;
end
save result222

% 求得滚子与套圈的pv值！！
for i=1:loadj
    pvzhi1(i)=Ph1ii(i)*deltaU1(i);
    pvzhi2(i)=Ph2ii(i)*deltaU2(i);
end
pvzhinonload=Ph1nonload*deltaU1nonload;
save T1;save T2;save T1nonload; save pvzhi1;save pvzhi2;save pvzhinonload;
save oilh1;save oilh2;save Fm;  save deltaU1nonload;save deltaU1;save deltaU2;save Wo;save Wx
