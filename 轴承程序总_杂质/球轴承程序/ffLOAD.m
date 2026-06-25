function ffLOAD(datafromvb,loadi,ttt3)  % 求解 内圈平衡 方程组

pp=ttt3;

% 数据输入部分！！
n=datafromvb(1);Dw=datafromvb(2);Dm=datafromvb(3);f1=datafromvb(4);f2=datafromvb(5);a0=datafromvb(6)*pi/180;
Rp=datafromvb(7); %兜孔半径
aqj=datafromvb(35)*pi/180; %安全角
ballden=datafromvb(61);e1=datafromvb(9);e2=datafromvb(10);e3=datafromvb(11);o1=datafromvb(12);o2=datafromvb(13);o3=datafromvb(14);
W2=datafromvb(15)*pi/30;Fxx=datafromvb(16); Fyy=datafromvb(17); Fzz=datafromvb(18); Myy=datafromvb(19)+0.1; Mzz=datafromvb(20)+0.1;
dangbianxishu1=datafromvb(21);dangbianxishu2=datafromvb(22);
cucao1=datafromvb(23);cucao2=datafromvb(24);cucao3=datafromvb(25);
%润滑油基本参数
oilden=datafromvb(26); sita0=datafromvb(27);K=datafromvb(28); %常温下的导热系数
niandu0=datafromvb(29); nianya0=datafromvb(30); beita0=datafromvb(31);   %粘温系数（近似认为不变，wys论文中）   所有参数均为常温下的参数！！！
yindao=datafromvb(32); yindaojianxi=datafromvb(33);   % 引导方式,引导间隙
dianpianjiao=datafromvb(34);  %垫片角

%a0=acos(1-deltar0/(f1+f2-1)/Dw/2);
deltac0=Rp/2-Dw/2;
deltar0=2*(1-cos(a0))*(f1+f2-1)*Dw;               %求原始径向间隙！
deltaa0=2*sin(a0)*(f1+f2-1)*Dw ;                 %求原始轴向间隙！
Dy1=Dm+Dw+deltar0;Dy2=Dm-Dw-deltar0;  % 沟底直径
Dr1=Dy1-2*dangbianxishu1*Dw;Dr2=Dy2+2*dangbianxishu2*Dw;  % 挡边直径（引导面直径）  档边高系数给出！！！！
save Dy1;
save Dy2;

a0=acos(1-deltar0/(2*(f1+f2-1)*Dw)-(2*f2-1)*(1-cos(dianpianjiao*pi/180))/(2*(f1+f2-1)));  %存在垫片角ad时计算初始接触角
xh=sqrt( (f2*Dw)^2-( f2*Dw*cos(dianpianjiao*pi/180)- dangbianxishu2*Dw )^2 );

z7=0;z8=0;z9=0;z10=0;z11=0;z12=0;z13=0;z14=0;z15=0;
loadj=length(loadi);
%初始的假设值！！！
beitajiao=atan(sin(a0)/(cos(a0)+Dw/Dm));
Wo=W2*(1-Dw/Dm*cos(a0))*(cos(a0)+tan(beitajiao)*sin(a0))/(( 1-Dw/Dm*cos(a0))*(cos(a0)+tan(beitajiao)*sin(a0))+ ( 1+Dw/Dm*cos(a0))*(cos(a0)+tan(beitajiao)*sin(a0))    );
Wx=-W2*(1-Dw/Dm*cos(a0))*(1+Dw/Dm*cos(a0))/(( 1-Dw/Dm*cos(a0))*(cos(a0)+tan(beitajiao)*sin(a0))+ ( 1+Dw/Dm*cos(a0))*(cos(a0)+tan(beitajiao)*sin(a0)) )/(Dw/Dm);

for i=1:loadj
    x(i)=pp(i);
    r(i)=pp(loadj+i)+Dm/2;
end
X2=pp(2*loadj+1);Y2=pp(2*loadj+2);Z2=pp(2*loadj+3);sitay=pp(2*loadj+4);sitaz=pp(2*loadj+5);


m=ballden*Dw*Dw*Dw*pi/6;  %球的质量计算
J=m*Dw*Dw/10;             %转动惯量计算
%求外、内圈赫兹接触载荷！
rou1=1/Dw*(4-1/f1-2*Dw*cos(a0)/(Dm+Dw*cos(a0)));
Rx1=Dw*(0.5*Dm/cos(a0)+0.5*Dw)/( Dm/cos(a0) );
Ry1=f1*Dw/(2*f1-1);
R21=Rx1*Ry1/(Rx1+Ry1);
E1=2/(((1-o1^2)/e1)+((1-o3^2)/e3));             %e1，e3为外圈和球的弹性模量，o1，o3为外圈，球的泊松比
K1=1.0339*(Ry1/Rx1)^0.636;
ee21=1.0003+0.5968/(Ry1/Rx1);          %ee21是外圈的ee2第二类完全椭圆积分p12 pengbo论文eq.2-16
ee11=((K1^2+1)*ee21-(K1^2-1)*ee21*( -2*cos(a0)/(Dm+Dw*cos(a0)) +1/(f1*Dw)  )/rou1)/2;

rou2= 1/Dw*(4-1/f2+2*Dw*cos(a0)/(Dm-Dw*cos(a0)));
E2=2/(((1-o2^2)/e2)+((1-o3^2)/e3));   %e1，e3为内圈和球的弹性模量，o1，o3为内圈，球的泊松比
Rx2=Dw*(0.5*Dm/cos(a0)-0.5*Dw)/(Dm/cos(a0) );
Ry2=f2*Dw/(2*f2-1);
K2=1.0339*(Ry2/Rx2)^0.636;
R22=Rx2*Ry2/(Rx2+Ry2);
R222=Dm/2+(f2-0.5)*Dw*cos(a0);
ee22=1.0003+0.5968/(Ry2/Rx2);  %ee21是内圈的ee2第二类完全椭圆积分p12 pengbo论文eq.2-16
ee12=((K2^2+1)*ee22-(K2^2-1)*ee22*( 2*cos(a0)/(Dm-Dw*cos(a0)) +1/(f2*Dw)  )/rou2)/2;
aaa1=0;                              %球和外圈实际接触角
aaa2=0;                              %球和内圈实际接触角

for i=1:loadj       %共loadj个球循环
    sita(i)=2*pi*(loadi(i))/n ;     %求各个球的方位角

    % 最小油膜厚度计算!!!
    if Fxx==0 || Fxx/((-Fzz*tan(a0)/Fxx+1)*n*sin(a0))<=0
        Qmax=5*Fzz/n;
    else
        Qmax=Fxx/((-Fzz*tan(a0)/Fxx+1)*n*sin(a0));
    end
    Q20base=max(0, 1-1/(2* (-2.5*Fzz*tan(a0)/max(abs(Fxx),eps)+2.5) )*(1-cos(sita(i))));
    Q20(i)=Qmax*Q20base^1.5;
    Q10(i)=Q20(i)+m*(Wo)^2*r(i);
    if Q10(i)<=0
        Q10(i)=1;
    end
    if Q20(i)<=0
        Q20(i)=1;
    end
    U1(i)=abs(-Wo*Dm/4-0.5*(abs(Wx)*cos(a0)+ Wo*cos(a0))*(Dw/2));
    U2(i)=abs((W2-Wo)*Dm/4+0.5*(abs(Wx)*cos(a0)-(W2-Wo)*cos(a0))*(Dw/2));    micro_config = load_micro_interface_config();
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
    oilh1(i)=thermal_factor1*texture_factor*Rx1*3.63*(niandu0*U1(i)/(E1*Rx1))^(0.68)*(nianya0*E1)^(0.49)*(Q10(i)/Rx1^2/E1)^(-0.073)*(1-exp(-0.68*K1));
    oilh2(i)=thermal_factor2*texture_factor*Rx2*3.63*(niandu0*U2(i)/(E2*Rx2))^(0.68)*(nianya0*E2)^(0.49)*(Q20(i)/Rx2^2/E2)^(-0.073)*(1-exp(-0.68*K2));

% oilh1(1)=0.114*10^(-6);oilh1(2)=0.114*10^(-6);oilh1(3)=0.114*10^(-6);oilh1(4)=0.114*10^(-6);oilh1(5)=0.115*10^(-6);oilh1(6)=0.116*10^(-6);
% oilh1(7)=0.117*10^(-6);oilh1(8)=0.118*10^(-6);oilh1(9)=0.118*10^(-6);
% oilh1(10)=0.119*10^(-6);oilh1(11)=0.119*10^(-6);oilh1(12)=0.119*10^(-6);oilh1(13)=0.118*10^(-6);oilh1(14)=0.118*10^(-6);oilh1(15)=0.117*10^(-6);
% oilh1(16)=0.117*10^(-6);oilh1(17)=0.116*10^(-6);oilh1(18)=0.115*10^(-6);
% 
% oilh2(1)=0.113*10^(-6);oilh2(2)=0.114*10^(-6);oilh2(3)=0.113*10^(-6);oilh2(4)=0.114*10^(-6);oilh2(5)=0.115*10^(-6);oilh2(6)=0.116*10^(-6);
% oilh2(7)=0.117*10^(-6);oilh2(8)=0.118*10^(-6);oilh2(9)=0.118*10^(-6);
% oilh2(10)=0.119*10^(-6);oilh2(11)=0.119*10^(-6);oilh2(12)=0.119*10^(-6);oilh2(13)=0.118*10^(-6);oilh2(14)=0.118*10^(-6);oilh2(15)=0.117*10^(-6);
% oilh2(16)=0.116*10^(-6);oilh2(17)=0.115*10^(-6);oilh2(18)=0.114*10^(-6);
    sitajiao(i)=360*(loadi(i))/n;
    P2(i)=5e-6*sin((Dw/2*sitajiao(i)/(0.001667))/180*pi+pi/2);%Chenzhegai 波纹度    % load-dependent debris_shift comes from micro_config when debris.enabled is true.
    xx(i)=x(i)+((f1-0.5)*Dw )*sin(a0);  y(i)=r(i)-(Dm/2-(f1-0.5)*Dw*cos(a0));%x(i)，r(i)变量
    %由位移-变形相容方程求弹性变形和接触角！
    %delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-0.542*10^(-6);
    %delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-0.542*10^(-6);
    delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift;%崔立博士论文eq.2-15，接触变形中涉及未知变量x(i)，r(i)，内圈接触还涉及X2，Y2,Z2,sitaz,sitay
    delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;

    % 如果 接触变形为负值，此时要改变滚动体位移！！

    if delta1(i)<=0
        xielv1=y(i)/xx(i);
        xielv=((f1+f2-1)*Dw*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i)))/( (f1+f2-1)*Dw*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i))));
        if xielv1>xielv
        else
            xxi=xx(i);
            xx(i)=2*( (y(i)*xielv+xx(i))/(xielv^2+1) )-xx(i);
            y(i)=2*xielv*( (y(i)*xielv+xxi)/(xielv^2+1) )-y(i);     %关于中心线作对称处理以保证斜率〉k2
            xielv1=y(i)/xx(i);
        end
    end
    while (delta1(i)<=0)%迭代 步长1e-7
        increase1=1e-7;
        xx(i)=xx(i)+increase1;
        x(i)=xx(i)-((f1-0.5)*Dw  )*sin(a0);
        y(i)=y(i)+increase1*abs(xielv1);
        r(i)=y(i)+(Dm/2-((f1-0.5)*Dw  )*cos(a0));
        delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw;
    end

   % delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-0.542*10^(-6);
   % delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-0.502*10^(-6);
     delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift;
     delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;

    if  delta2(i)<=0
        xielv=((f1+f2-1)*Dw*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i)))/( (f1+f2-1)*Dw*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i))));
        if (y(i)/xx(i))>xielv
        else
            xxi=xx(i);
            xx(i)=2*( (y(i)*xielv+xx(i))/(xielv^2+1) )-xx(i);
            y(i)=2*xielv*( (y(i)*xielv+xxi)/(xielv^2+1) )-y(i);     %关于中心线作对称处理以保证斜率〉k2
        end
        xielv2=(  (((f1+f2-1)*Dw  )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))) -y(i) )/(   (((f1+f2-1)*Dw  )*sin(a0)+(X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i))))) -xx(i));
    end

    %delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-0.4*10^(-6);
    %delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-0.502*10^(-6);
    delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift;
    delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;

    while (delta1(i)<=0)
        increase1=1e-7;
        xx(i)=xx(i)+increase1;%调节 r(i)
        x(i)=xx(i)-((f1-0.5)*Dw  )*sin(a0);
        y(i)=y(i)+increase1*xielv1;%调节 r(i)
        r(i)=y(i)+(Dm/2-((f1-0.5)*Dw  )*cos(a0));
        delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw -oilh1(i)+debris_shift;
    end

    while (delta2(i)<=0)
        if  xx(i)<1e-10   |   sqrt(xx(i)^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift<1e-10                    %如果xx(i)<0, 要停止改变x(i) r(i),而是改变X2
            while(delta2(i)<=0)
                X2=X2+1e-7;%调节X2
                delta2(i)=sqrt(( ((f1+f2-1)*Dw    )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw    )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;
            end
        end
        decrease1=1e-7;
        xx(i)=xx(i)-1e-7;
        x(i)=xx(i)-((f1-0.5)*Dw  )*sin(a0);
        y(i)=y(i)-decrease1*xielv2;
        r(i)=y(i)+(Dm/2-((f1-0.5)*Dw  )*cos(a0));
        delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;

        if   sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift<1e-10            %如果改变xx(i) y(i)使得delta1(i)<0 则要停止改变x(i) r(i),而是改变X2
            delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift;
            if delta1(i)<=0
                xielv1=y(i)/xx(i);
                xielv=((f1+f2-1)*Dw*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i)))/( (f1+f2-1)*Dw*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i))));
                if xielv1>xielv
                else
                    xxi=xx(i);
                    xx(i)=2*( (y(i)*xielv+xx(i))/(xielv^2+1) )-xx(i);
                    y(i)=2*xielv*( (y(i)*xielv+xxi)/(xielv^2+1) )-y(i);     %关于中心线作对称处理以保证斜率〉k2
                    xielv1=y(i)/xx(i);
                end
            end

            while (delta1(i)<=0)
                increase1=1e-7;
                xx(i)=xx(i)+increase1;
                x(i)=xx(i)-((f1-0.5)*Dw  )*sin(a0);
                y(i)=y(i)+increase1*xielv1;
                r(i)=y(i)+(Dm/2-((f1-0.5)*Dw  )*cos(a0));
                delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw -oilh1(i)+debris_shift;
            end
            delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;
            while(delta2(i)<=0)
                X2=X2+1e-7;
                delta2(i)=sqrt(( ((f1+f2-1)*Dw    )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw    )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw ) -oilh2(i)+debris_shift;
                delta2(i)=sqrt(( ((f1+f2-1)*Dw    )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw    )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw ) -oilh2(i)+debris_shift;

            end
        end
    end

    %delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-0.4*10^(-6);
    %delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-0.542*10^(-6);
     delta1(i)=sqrt((xx(i))^2+y(i)^2)-(f1-0.5)*Dw-oilh1(i)+debris_shift;
     delta2(i)=sqrt(( ((f1+f2-1)*Dw )*sin(a0)+X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i)))-(xx(i)))^2+( ((f1+f2-1)*Dw )*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i))-y(i))^2)-((f2-0.5)*Dw )-oilh2(i)+debris_shift;

    Q1(i)=sqrt( ( delta1(i)*pi*(2*ee21*K1^2/(pi))^0.333/(ee11))^3*4*E1^2/(9*rou1) );
    Q2(i)=sqrt( ( delta2(i)*pi*(2*ee22*K2^2/(pi))^0.333/(ee12))^3*4*E2^2/(9*rou2) );
    a1(i)=atan((xx(i))/y(i));
    a2(i)=atan((((f1+f2-1)*Dw    )*sin(a0)+(X2+R222*(sitaz*sin(sita(i))+sitay*cos(sita(i))))-(xx(i)))/(((f1+f2-1)*Dw)*cos(a0)+Z2*cos(sita(i))+Y2*sin(sita(i)) -y(i)));
    aa1(i)=(6*K1^2*ee21*R21*Q1(i)/(E1*pi))^0.3333; b1(i)=aa1(i)/K1;    %求得接触椭圆长短轴！
    aa2(i)=(6*K2^2*ee22*R22*Q2(i)/(E2*pi))^0.3333; b2(i)=aa2(i)/K2;   %求得接触椭圆长短轴！

    jcqwj(i)=aa1(i)/Dw;           %每个球外圈接触椭圆与滚动体质心的最大张角
    jcqnj(i)=aa2(i)/Dw;           %每个球内圈接触椭圆与滚动体质心的最大张角

    if a1(i)+jcqwj(i)>aaa1               %求球和外圈实际接触角
        aaa1=a1(i)+jcqwj(i);
    end
    if a2(i)+jcqnj(i)>aaa2               %求球和内圈实际接触角
        aaa2=a2(i)+jcqnj(i);
    end

    Ph1(i)=1.5*Q1(i)/(pi*aa1(i)*b1(i));
    Ph2(i)=1.5*Q2(i)/(pi*aa2(i)*b2(i));

    %求球沿椭圆长轴方向的摩擦力！ 最初取摩擦系数为0.02，待转速求得后再重新收敛！！！！
    deltaU1(i)=abs(-Wo*Dm/2+abs(Wx)*cos(a0)-Wo*cos(a0)*(Dw/2));
    S12(i)=abs(deltaU1(i)/U1(i));  
    miuI(i)=0.0127*(50/(50-S12(i)))*log(0.584*Q1(i)/niandu0/deltaU1(i)/(U1(i))^2);    %CHENZHEGAI
    deltaU2(i)=abs(+(W2-Wo)*Dm/2- abs(Wx)*cos(a0)+(W2-Wo)*cos(a0)*(Dw/2));
    S22(i)=abs(deltaU2(i)/U2(i));
    miuO(i)=0.0127*(50/(50-S22(i)))*log(0.584*Q2(i)/niandu0/deltaU2(i)/(U2(i))^2);    %CHENZHEGAI
    fs1(i)=miuI(i)*Q1(i);fs2(i)=miuO(i)*Q2(i);

    Fz(i)=m*(Wo)^2*r(i);
    z1(i)=-Q1(i)*sin(a1(i))+Q2(i)*sin(a2(i))+fs1(i)*cos(a1(i))-fs2(i)*cos(a2(i));%球x方向受力平衡(a)，fs1与崔力论文中计算不同
    z3(i)=-Q1(i)*cos(a1(i))+Q2(i)*cos(a2(i))-fs1(i)*sin(a1(i))+fs2(i)*sin(a2(i))+Fz(i);%球x方向受力平衡(c)

    if yindao==2   %内圈引导
        %求引导面/保持架之间的法向力和摩擦力！
        Kbp=11/yindaojianxi;
        deltayindao=sqrt(Y2^2+Z2^2);     % 偏离原来位置的距离！
        Fyindao=Kbp*abs(deltayindao);
        fyindao=0.05*Fyindao;
        zuoyongjiao=atan(Y2/Z2);    save Fyindao;save fyindao;save zuoyongjiao
        z12=z12+(Q2(i)*cos(a2(i))+fs2(i)*sin(a2(i)))*sin(sita(i))+Fyindao*sin(zuoyongjiao);%内圈平衡方程（c）
        z13=z13+(Q2(i)*cos(a2(i))+fs2(i)*sin(a2(i)))*cos(sita(i))+Fyindao*cos(zuoyongjiao);%内圈平衡方程（d）
    else     %外圈引导 引导面对内圈没有作用力
        z12=z12+(Q2(i)*cos(a2(i))+fs2(i)*sin(a2(i)))*sin(sita(i));
        z13=z13+(Q2(i)*cos(a2(i))+fs2(i)*sin(a2(i)))*cos(sita(i));
        Fyindao=0;    fyindao=0;        zuoyongjiao=0;
    end
    save Fyindao;save fyindao;save zuoyongjiao

    z11=z11+(Q2(i)*sin(a2(i))-fs2(i)*cos(a2(i)));%内圈平衡方程（a）
    z14=z14+((Dm/2+(f2-0.5)*Dw*cos(a0))*(Q2(i)*sin(a2(i))-fs2(i)*cos(a2(i)))+f2*Dw*fs2(i)*cos(a2(i)) )*cos(sita(i));%内圈平衡方程（d）
    z15=z15+((Dm/2+(f2-0.5)*Dw*cos(a0))*(Q2(i)*sin(a2(i))-fs2(i)*cos(a2(i)))+f2*Dw*fs2(i)*cos(a2(i)) )*sin(sita(i));%内圈平衡方程（e）
end
Doc=Dy1-Dw*(1-cos(aaa1+aqj*pi/180));      %外滚道内表面直径
Dic=Dy2+Dw*(1-cos(aaa2+aqj*pi/180));      %外滚道内表面直径
save Doc;
save Dic;

z11=Fxx-z11;    z12=Fyy-z12;    z13=Fzz-z13;    z14=Myy-z14;    z15=Mzz-z15;

zzzz3=[z1';z3'; z11 ;z12;  z13; z14; z15 ];save zzzz3;            % 计算轴承 载荷用
result333=0;
for i=1:2*loadj+5
    result333=result333+zzzz3(i)^2
end
save result333;

result444=0;
for i=1:2*loadj+3
    result444=result444+zzzz3(i)^2;
end
save result444;

for i=1:loadj
    r(i)=r(i)-Dm/2;
end
uu=[ x r X2 Y2 Z2 sitay sitaz ]';            %当位移变量有所改变时，需要对变量更新！！
save uu

q1q2a1a2=[Q1 Q2 a1 a2];
save q1q2a1a2
save Ph1; save Ph2; save aa1;save aa2;save Fz;

% 计算接触点与内圈挡边的夹角！！！
for i=1:loadj
    jiaodu2(i)=atan( (xh-X2+x(i))/(Dm/2-Dr2/2+r(i)))-a2(i);
end
jiaodu2=jiaodu2*180/pi;
save jiaodu2

