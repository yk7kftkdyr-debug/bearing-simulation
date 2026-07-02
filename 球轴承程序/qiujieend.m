function [loadj,returndata]=qiujieend(datafromvb,micro_config)

%定义如下参数
%1-滚动体数；2-滚动体直径；3-节圆直径；4-外圈曲率系数；5-内圈曲率系数；6-接触角；7-兜孔半径；35-安全角；
%61-滚动体密度；68-保持架密度；51-外圈密度；50-内圈密度；9-内圈弹性模量；10-外圈弹性模量；11-滚动体弹性模量；
%65-保持架弹性模量；12-内圈泊松比；13-外圈泊松比；14-球泊松比；66-保持架泊松比；69-外圈转速；15-内圈转速；16-轴载；17-Y载
%18-Z载荷;19-Y弯矩；20-Z弯矩；21-外圈挡边系数；22-内圈挡边系数；23-内圈粗糙度；24-外圈粗糙度；25-球粗度度；26-油密度
%27-入口油温；28-导热系数；29-粘度；30-粘压系数；31-粘温系数；32-引导方式；33-引导间隙；34-垫片角；
%74-轴承类型；36-外圈直径；37-内圈直径；38-轴直径；39-轴承座直径；40-轴弹性模量；41-轴承座弹模；42-轴泊松比；43-轴承座泊松比    
%44-轴承/轴承座过盈量；45-轴/轴承过盈；54-轴线膨胀系数；55-轴承座线膨胀；56-内圈线膨胀；57-外圈线膨胀；58-球线膨胀；46-球温
%47-内圈温度；48-外圈温度；49-环境温度；59-轴温度；60-轴承座温度；50-内圈密度；51-外圈密度；52-轴密度；53-轴承座密度
%70-轴承润滑方式；71-球/套圈摩擦系数；72-球/保持架摩擦系数；73-保持架/套圈摩擦系数；  

if nargin < 1 || isempty(datafromvb)
datafromvb=zeros(75,1);
datafromvb(1)=15;datafromvb(2)=0.02223;datafromvb(3)=0.1253;datafromvb(4)=0.5232;datafromvb(5)=0.5232;datafromvb(6)=40; datafromvb(7)=0.0115;datafromvb(35)=5;
datafromvb(61)=7800; datafromvb(68)=7800; datafromvb(51)=7800; datafromvb(50)=7800; datafromvb(9)=2.06e+11;datafromvb(10)=2.06e+11;datafromvb(11)= 2.06e+11; 
datafromvb(65)=2.06e+11; datafromvb(12)=0.3; datafromvb(13)=0.3; datafromvb(14)=0.3; datafromvb(66)=0.3;datafromvb(69)=0; datafromvb(15)=10000;datafromvb(16)=20000 ;datafromvb(17)=0; 
datafromvb(18)=3000; datafromvb(19)=0; datafromvb(20)=0;datafromvb(21)=0.275;datafromvb(22)=0.225;datafromvb(23)=1e-07;datafromvb(24)=1e-07;datafromvb(25)=1e-07;datafromvb(26)=970; 
datafromvb(27)=20;datafromvb(28)=0.0966;datafromvb(29)=0.0318; datafromvb(30)=1.28e-008; datafromvb(31)=3.2e-002; datafromvb(32)=1; datafromvb(33)=0.4e-03; datafromvb(34)=3.8; 
datafromvb(74)=2;datafromvb(36)=370e-03;datafromvb(37)=220e-03;datafromvb(38)=120.0e-03;datafromvb(39)=225e-03;datafromvb(40)=1.96e11;datafromvb(41)=2.18e11;datafromvb(42)=0.3;datafromvb(43)=0.3;
datafromvb(44)=0.5e-03;datafromvb(45)=-0.5e-03;datafromvb(54)=11.6e-06;datafromvb(55)=11.8e-06;datafromvb(56)=11.8e-06;datafromvb(57)=11.8e-06;datafromvb(58)=11.8e-06;datafromvb(46)=180;
datafromvb(47)=170;datafromvb(48)=190;datafromvb(49)=27;datafromvb(59)=195;datafromvb(60)=160;datafromvb(50)=7870;datafromvb(51)=7870;datafromvb(52)=7860;datafromvb(53)=8360;
datafromvb(70)=0; datafromvb(71)=0.2; datafromvb(72)=0.2; datafromvb(73)=0.2;
end

projectRoot=fileparts(fileparts(mfilename('fullpath')));
if exist(projectRoot,'dir')==7
    addpath(projectRoot);
end
if nargin < 2 || isempty(micro_config)
    micro_config = make_micro_interface_config();
end
save('micro_config_runtime.mat','micro_config');

data_input=datafromvb;
[deltaw_input, clearInfo] = calcBallWorkingClearance(data_input);
Dw_input=data_input(2); f1_input=data_input(4); f2_input=data_input(5);
cos_a_work=1-deltaw_input/(2*(f1_input+f2_input-1)*Dw_input);
if cos_a_work>1
    cos_a_work=1;
elseif cos_a_work<-1
    cos_a_work=-1;
end
data_work=data_input;
data_work(6)=acos(cos_a_work)*180/pi;
datafromvb=data_work;
fprintf('\n===== 球轴承游隙检查 =====\n');
fprintf('真实输入初始游隙 = %.6e m = %.6f mm\n', 2*(1-cos(data_input(6)*pi/180))*(data_input(4)+data_input(5)-1)*data_input(2), 2*(1-cos(data_input(6)*pi/180))*(data_input(4)+data_input(5)-1)*data_input(2)*1000);
fprintf('参与求解工作游隙 = %.6e m = %.6f mm\n', deltaw_input, deltaw_input*1000);
fprintf('==========================\n');

% 数据输入部分！！
n=datafromvb(1);Dw=datafromvb(2);Dm=datafromvb(3);f1=datafromvb(4);f2=datafromvb(5);a0=datafromvb(6)*pi/180;
Rp=datafromvb(7); %兜孔半径
aqj=datafromvb(35)*pi/180; %安全角
ballden=datafromvb(61);e1=datafromvb(9);e2=datafromvb(10);e3=datafromvb(11);o1=datafromvb(12);o2=datafromvb(13);o3=datafromvb(14); 
W2=datafromvb(15)*2*pi/60;Fxx=datafromvb(16); Fyy=datafromvb(17); Fzz=datafromvb(18); Myy=datafromvb(19)+0.1; Mzz=datafromvb(20)+0.1; 
dangbianxishu1=datafromvb(21);dangbianxishu2=datafromvb(22);
cucao1=datafromvb(23);cucao2=datafromvb(24);cucao3=datafromvb(25); 
%润滑油基本参数
oilden=datafromvb(26); sita0=datafromvb(27);K=datafromvb(28); %常温下的导热系数
niandu0=datafromvb(29); nianya0=datafromvb(30); beita0=datafromvb(31);   %粘温系数（近似认为不变，wys论文中）   所有参数均为常温下的参数！！！
yindao=datafromvb(32); yindaojianxi=datafromvb(33);   % 引导方式,引导间隙
dianpianjiao=datafromvb(34)*pi/180;  %垫片角
leixing=datafromvb(74); %轴承类型
deltar0=2*(1-cos(a0))*(f1+f2-1)*Dw;               %求原始径向间隙！
%a0=acos(1-deltar0/(f1+f2-1)/Dw/2);
qiujieLOAD(datafromvb);  
load loadi;load www3          %初步求解载荷分布！！
ffLOAD(datafromvb,loadi,www3);load zzzz3;zz3=zzzz3;
iii=1;j=1;
for iii=1:5
    neiquan1=www3; neiquan1(2*loadj+iii)=www3(2*loadj+iii)+1e-10;
    ffLOAD(datafromvb,loadi,neiquan1);load zzzz3;
    for j=1:5
        kk(iii,j)=(zzzz3(2*loadj+j)-zz3(2*loadj+j))/1e-10;
    end
end
disp('刚度矩阵为：')
kk=-kk';
save kk
ffLOAD(datafromvb,loadi,www3);load q1q2a1a2;

qiujieSPEED1(datafromvb); 
load wwmin34;
for i=1:loadj
   Wo(i)=wwmin34(i);      Wx(i)=wwmin34(loadj+i);    Wz(i)=wwmin34(2*loadj+i);  Wy(i)=-Wo(i)/100;
end
if loadj==n
    Wononload=0; Wxnonload=0; Wznonload=0; Wynonload=0;
else
    Wononload=wwmin34(3*loadj+1); Wxnonload=wwmin34(3*loadj+2); Wznonload=0; Wynonload=0;
end

for i=1:loadj
beitajiao(i)=atan(sin(a1(i))/(cos(a1(i))+Dw/Dm));
Wolilun(i)=W2*(1-Dw/Dm*cos(a2(i)))*(cos(a1(i))+tan(beitajiao(i))*sin(a1(i)))/(( 1-Dw/Dm*cos(a2(i)))*(cos(a1(i))+tan(beitajiao(i))*sin(a1(i)))+ ( 1+Dw/Dm*cos(a1(i)))*(cos(a2(i))+tan(beitajiao(i))*sin(a2(i)))    );
end
lilunWc=max(Wolilun);
zhuansuxishu=1;
if  lilunWc<max(Wo)
    zhuansuxishu=lilunWc/max(Wo);
end
Wo=Wo*zhuansuxishu;     Wononload=Wononload*zhuansuxishu;
disp('保持架转速为：理论转速：打滑率为：')
Wc=(sum(Wo)+ Wononload*(n-loadj))/n; 
if Wc>lilunWc
    dahualv=1e-3;           %此时认为不打滑，！！
else
    dahualv=(1-Wc/lilunWc);
end

load pvzhi1;load pvzhi2;load pvzhi1nonload;
%旋滚比计算！！
for i=1:loadj
    sitajiao(i)=360*(loadi(i))/n; 
    beita(i)=atan(sin(a1(i))/(cos(a1(i))+Dw/Dm  ));
    xuan1(i)= -(1+Dw/Dm*cos(a1(i)))*tan(a1(i)-beita(i))+Dw/Dm*sin(a1(i))   ;
    xuan2(i)= (1-Dw/Dm*cos(a2(i)))*tan(a2(i)-beita(i))+Dw/Dm*sin(a2(i))   ;
end
save xuan1;xuan2;
allsitajiao=360/n:360/n:360;  m=ballden*Dw*Dw*Dw*pi/6;  %球的质量计算
for i=1:n
    for j=1:loadj
      if i==loadi(j)
          mark1=i;break
      else 
          mark1=0;
      end
    end
    if mark1==i
        allQ1(i)=Q1(j);
        allQ2(i)=Q2(j);
        allPh1(i)=Ph1(j);
        allPh2(i)=Ph2(j);
        allpvzhi1(i)=pvzhi1(j);
        allpvzhi2(i)=pvzhi2(j);
        allWo(i)=Wo(j);
        allWx(i)=Wx(j);
        allxuan2(i)=xuan2(j);
    else
        allQ1(i)=m*Wononload^2*Dm/2;
        allQ2(i)=0; allPh2(i)=0;
        aa11=(6*K1^2*ee21*R21*allQ1(i)/(E1*pi))^0.3333; b11=aa11/K1;    %求得接触椭圆长短轴！
        allPh1(i)=1.5*allQ1(i)/(pi*aa11*b11);
        allpvzhi1(i)=pvzhi1nonload;
        allpvzhi2(i)=0;
        allWo(i)=Wononload;
        allWx(i)=Wxnonload;
        allxuan2(i)=0;
    end
end

Famax1(data_input);load Famax;load Fyindao;load fyindao;load zuoyongjiao
deltapd=clearInfo.deltapd; deltapt1=clearInfo.deltapt1; deltapt2=clearInfo.deltapt2; deltapf1=clearInfo.deltapf1; deltapf2=clearInfo.deltapf2; deltaw=clearInfo.deltaw;
load deltaU1;load deltaU2;load T1;load T2;load jiaodu2;
load oilh1;load oilh2;                 %最小油膜厚度！！

life_ball(datafromvb,Q2,Q1,a2,a1,leixing,Wx)
load life_ball Lgundao L1 L2 L3

%非承载半圈最小间隙
for i=1:loadj
    orv=sqrt(((2*f2-1)*Dw*sin(dianpianjiao))^2+(Dm/2+(f2-0.5)*Dw-deltar0/4)^2);
    racecenterx=((f2-0.5)*Dw*sin(dianpianjiao)-www3(2*loadj+1)+orv*sin(www3(2*loadj+4))*cos(sitajiao(i))-orv*sin(www3(2*loadj+5))*sin(sitajiao(i))-www3(i));
    racecentery=Dm/2+(f2-0.5)*Dw-deltar0/4+sqrt((www3(2*loadj+2)+orv*(cos(www3(2*loadj+5))-1)*sin(sitajiao(i)))^2*(sin(sitajiao(i)))^2+(www3(2*loadj+3)+orv*(cos(www3(2*loadj+4))-1)*cos(sitajiao(i)))^2*(cos(sitajiao(i)))^2)-(www3(loadj+i)+Dm/2);
    beitaa=atan(racecentery/racecenterx);
    dist=sqrt(racecentery^2+racecenterx^2+Dw^2/4-sqrt(racecentery^2+racecenterx^2)*Dw*cos(3/2*pi-beitaa-a2(i)));
    dmin(i)=f2*Dw-dist;
end

hminmin=min(oilh2);

for i=1:loadj   
    beita1(i)=acos(1-2*(2*aa1(i))^2/Dw^2);    
    dangbianxishu11(i)=1-cos(a1(i)+beita1(i)/2+aqj);    
    Drr1(i)=Dy1-2*dangbianxishu11(i)*Dw;
    
    beita2(i)=acos(1-2*(2*aa2(i))^2/Dw^2);    
    dangbianxishu22(i)=1-cos(a2(i)+beita2(i)/2+aqj);    
    Drr2(i)=Dy2+2*dangbianxishu22(i)*Dw;  
end

Dr1new=min(Drr1); Dr2new=max(Drr2);
dangbianxishu1new=(Dy1-Dr1new)/2/Dw;    dangbianxishu2new=(Dr2new-Dy2)/2/Dw;
for i=1:loadj
    jxj1(i)=acos(1-dangbianxishu1)-a1(i)-beita1(i)/2;
    jxj2(i)=acos(1-dangbianxishu2)-a2(i)-beita2(i)/2;
    aqj1(i)=acos(1-dangbianxishu1new)-a1(i)-beita1(i)/2;
    aqj2(i)=acos(1-dangbianxishu2new)-a2(i)-beita2(i)/2;
end



% 写入文本用来保存
% 页眉
fid=fopen(bearingOutputPath('bearing','qiu','ballbearing.txt'),'w','n','UTF-8');
fprintf(fid,'\t\t\t\t\t\t************************\n');
fprintf(fid,'\t\t\t\t\t\t高速滚动轴承拟动力学分析 \n');
fprintf(fid,'\t\t\t\t\t\t************************\n');
% 输入部分
fprintf(fid,'\n*******************************输 入 参 数**********************************\n');
fprintf(fid,'\n 轴承中径（mm）  滚动体直径（mm）   滚动体数目     接触角(°)      内圈曲率系数    外圈曲率系数     安全角(°)     \n');
fprintf(fid,'%9d   %15.3d  %11d  %12.3d  %15.3d  %13.3d  %13.3d \n', [ Dm*1000    Dw*1000   n     a0*180/pi    f2   f1   aqj]);
fprintf(fid,'\n************************************************************************');
fprintf(fid,'\n            弹 性 模 量 (Pa)                     泊松比   \n');
fprintf(fid,'   外圈         内圈          球          外圈         内圈          球 \n');
fprintf(fid,'%11.3d   %11.3d  %11.3d  %11.3d  %11.3d %11.3d  %11.3d  \n', [ e1  e2   e3   o1  o2  o3 ]);
fprintf(fid,'\n\n************************************************************************');
fprintf(fid,'\n 滑油密度（kg/m3）  入口油温（℃）   热传导系数     粘度(Pa*s)        粘压系数     粘温系数  \n');
fprintf(fid,'%9d   %15.3d  %17.3d  %12.3d  %11.3d  %11.3d  \n', [ oilden    sita0  K  niandu0   nianya0  beita0 ]);
fprintf(fid,'\n************************************************************************');
fprintf(fid,'\n  内圈转速                    轴上的负荷   \n');
fprintf(fid,'  r/min            沿 X(N)       沿 Y(N)      沿 Z(N)       绕 Y(N*m)       绕 Z(N*m)   \n');
fprintf(fid,'%7d   %9d  %9d  %9d  %9d  %9d  \n', [ W2*30/pi  Fxx Fyy Fzz Myy-0.1 Mzz-0.1 ]);
% 计 算 结 果
fprintf(fid,'\n\n\n*****************************计 算 结 果***********************************\n');
fprintf(fid,'\n  球     球             接触负荷               接触角                 接触区域全长            接触区域全宽              沟底直径              挡边直径               最小油膜厚度\n');
fprintf(fid,'  号   方位角      外圈        内圈        外圈        内圈        外圈        内圈        外圈        内圈          外圈        内圈        外圈       内圈       外圈        内圈\n');
fprintf(fid,'       deg          N           N          deg         deg          m           m           m           m           mm         mm           mm         mm \n');
for i=1:loadj
if loadi(i)<10
    fprintf(fid,'\n %2d  %7.1f  %11.3d  %9.3d  %9.3d  %9.3d  %9.3d   %9.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d   %9.3d  %9.3d  %9.3d \n', [loadi(i) sitajiao(i)  Q1(i) Q2(i) a1(i)*180/pi a2(i)*180/pi aa1(i)*2 aa2(i)*2 b1(i)*2 b2(i)*2  Dy1   Dy2  Doc  Dic  oilh1(i)   oilh2(i)]);
else
    fprintf(fid,'\n %1d  %7.1f  %11.3d  %9.3d  %9.3d  %9.3d  %9.3d   %9.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d   %9.3d  %9.3d  %9.3d \n', [loadi(i) sitajiao(i)  Q1(i) Q2(i) a1(i)*180/pi a2(i)*180/pi aa1(i)*2 aa2(i)*2 b1(i)*2 b2(i)*2  Dy1   Dy2  Doc  Dic  oilh1(i)   oilh2(i)]);
end    
end
fprintf(fid,'\n\n*************************************************************************');
fprintf(fid,'\n  球     接触区椭圆对应中心角          原挡边极限角度        按安全角设计后的挡边直径        实际安全角 \n');
fprintf(fid,'  号        外圈        内圈           外圈        内圈         外圈        内圈         外圈        内圈     \n');
fprintf(fid,'            deg         deg           deg         deg           mm          mm         deg          deg     \n');
for i=1:loadj
if loadi(i)<10
    fprintf(fid,'\n %2d  %11.3f  %11.3f  %11.3f  %11.3f  %15.3f  %11.3f   %11.3f  %11.3f  \n', [loadi(i) beita1(i)*180/pi  beita2(i)*180/pi jxj1(i)*180/pi jxj2(i)*180/pi Dr1new*1000 Dr2new*1000 aqj1(i)*180/pi aqj2(i)*180/pi    ]);
else
    fprintf(fid,'\n %1d  %11.3f  %11.3f  %11.3f  %11.3f  %15.3f  %11.3f   %11.3f  %11.3f  \n', [loadi(i) beita1(i)*180/pi  beita2(i)*180/pi jxj1(i)*180/pi jxj2(i)*180/pi Dr1new*1000 Dr2new*1000 aqj1(i)*180/pi aqj2(i)*180/pi    ]);
end    
end
fprintf(fid,'\n\n*************************************************************************');
fprintf(fid,'\n     最大压应力               接触变形               自转          自转         自转        公转        离心力     非承载半圈最小间隙\n');
fprintf(fid,  '  外圈         内圈       外圈        内圈          角速度        角速度       角速度      角速度\n');
fprintf(fid,  '   Pa          Pa           m           m         绕X（r/min）  绕Y（r/min）  绕Z（r/min）  r/min       N          mm\n');
for i=1:loadj
    fprintf(fid,' %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d  %9.3d   %9.3d   %9.3d\n', [Ph1(i) Ph2(i) delta1(i) delta2(i) Wx(i)*30/pi  Wy(i)*30/pi  Wz(i)*30/pi  Wo(i)*30/pi    m*Wo(i)^2*Dm/2    dmin(i)*1000]);
end     
fprintf(fid,'\n\n*************************************************************************');
fprintf(fid,'\n  接触区最大滑动速度            PV值                    拖动力            旋滚比       球/兜孔     球位移    球位移\n');
fprintf(fid,  '  外圈         内圈       外圈        内圈        外圈        内圈        内圈         冲击力       轴向      径向   \n');
fprintf(fid,  '   m/s         m/s       Pam/s       Pam/s        N           N           m            N           m         m   \n');
for i=1:loadj
    fprintf(fid,' %1.3d  %9.3d  %9.3d  %9.3d  %9.3d   %9.3d  %9.3d %9.3d %9.3d %9.3d \n', [deltaU1(i) deltaU2(i) pvzhi1(i) pvzhi2(i) T1(i)  T2(i)  xuan2(i) Fsrx(i) www3(i) www3(loadj+i) ]);
end 
fprintf(fid,'\n\n*************************************************************************');
fprintf(fid,'\n                      内圈相对外圈的总位移 \n');
fprintf(fid,'    沿 X         沿 Y         沿 Z         绕 Y         绕 Z   \n');
fprintf(fid,'     m            m            m           min          min   \n');
fprintf(fid,'  %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  \n', [ www3(2*loadj+1)  www3(2*loadj+2) www3(2*loadj+3) www3(2*loadj+4)*180*60/pi  www3(2*loadj+5)*180*60/pi  ]);
fprintf(fid,'\n\n*************************************************************************');
fprintf(fid,'\n                刚度矩阵(作用力对位移的偏导数)\n');
fprintf(fid,'        x           y           z          deltay       deltaz\n');
fprintf(fid,'       N/m         N/m         N/m         N/rad        N/rad \n');
fprintf(fid,' Fx  %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  \n', [ kk(1,1)  kk(1,2) kk(1,3) kk(1,4) kk(1,5)  ]);
fprintf(fid,' Fy  %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  \n', [ kk(2,1)  kk(2,2) kk(2,3) kk(2,4) kk(2,5)  ]);
fprintf(fid,' Fz  %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  \n', [ kk(3,1)  kk(3,2) kk(3,3) kk(3,4) kk(3,5)  ]);
fprintf(fid,' My  %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  \n', [ kk(4,1)  kk(4,2) kk(4,3) kk(4,4) kk(4,5)  ]);
fprintf(fid,' Mz  %1.3d  %9.3d  %9.3d  %9.3d  %9.3d  \n', [ kk(5,1)  kk(5,2) kk(5,3) kk(5,4) kk(5,5)  ]);
fprintf(fid,'\n\n*************************************************************************');
if yindao==1
fprintf(fid,'\n   保持架打滑度      档边极限载荷     最小档边安全角   \n');
fprintf(fid,'  %8.1f%%    %14.3f    %14.3f   \n', [ dahualv*100   Famax  min(jiaodu2)  ]);
else
fprintf(fid,'\n   保持架打滑度      档边极限载荷     最小档边安全角      引导面法向力     引导面摩擦力  \n');
fprintf(fid,'   %8.1f%%      %14.3d     %14.3f     %9.3f    %9.3f  \n', [ dahualv*100   Famax  min(jiaodu2)  Fyindao  fyindao ]);   
end
fprintf(fid,'\n\n*************************************************************************');
fprintf(fid,'\n   装配引起的游隙变化(m)      温差引起的游隙变化(m)     离心力引起的游隙变化(m)      温差对配合影响引起的游隙变化(m)     离心力对配合影响引起的游隙变化(m)     工作游隙(m)\n');
fprintf(fid,'  %8.1d    %25.3d    %25.3d    %25.3d   %25.3d    %25.3d   \n', [ deltapd   deltapt1  deltapf1    deltapt2  deltapf2   deltaw]);
fclose(fid);;


%返回值
returndata=[kk(1,1)  kk(2,2)  kk(3,3)  kk(4,4)  kk(5,5)  dahualv  Doc  Dic  L3*1e6/(W2*30/pi*60) hminmin]; 


% 绘图供vb调用！！
set(gcf,'visible','off')
k=0.45; 
set(gcf,'units',get(gcf,'paperunits')); 
set(gcf,'paperposition',get(gcf,'position')*k);   % 改变 图片的尺寸而不影响质量   
set(gcf,'DefaultAxesFontSize',6)               % 设置 图片中的默认字体               
subplot(2,1,1); plot(allsitajiao,allQ1);xlabel('方位角（度）');ylabel('接触载荷（N）');subplot(2,1,2); plot(allsitajiao,allQ2);xlabel('方位角（度）');ylabel('接触载荷（N）');
print(gcf,'-djpeg',bearingOutputPath('bearing','qiu','loadball.jpg'))    % 球接触载荷分布
subplot(2,1,1); plot(allsitajiao,allpvzhi1);xlabel('方位角（度）');ylabel('PV值（Pam/s）');subplot(2,1,2); plot(allsitajiao,allpvzhi2);xlabel('方位角（度）');ylabel('PV值（Pam/s）');
print(gcf,'-djpeg',bearingOutputPath('bearing','qiu','pvball.jpg'))      % pv值分布
subplot(1,1,1);plot(allsitajiao,allxuan2);xlabel('方位角（度）');ylabel('旋滚比');
print(gcf,'-djpeg',bearingOutputPath('bearing','qiu','spintorollball.jpg'))      % 旋滚值分布
subplot(1,1,1);plot(allsitajiao,allWo*30/pi);xlabel('方位角（度）');ylabel('转速(r/min)');
print(gcf,'-djpeg',bearingOutputPath('bearing','qiu','rotatespeedball.jpg'))      % 公转转速分布

