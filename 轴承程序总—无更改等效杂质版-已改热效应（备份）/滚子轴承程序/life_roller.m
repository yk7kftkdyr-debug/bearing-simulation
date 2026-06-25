function life_roller(datafromvb,ROLLER1,ROLLER2,leixing,w1)
% ROLLER1 是内滚道与滚子的接触载荷分布
% ROLLER2 是外滚道与滚子的接触载荷分布
% aa是接触角（°）
% l是滚子“有效”长度
% D 是滚子直径
% Dm 是节圆直径
% leixing 是滚子轴承类型代号
% w1 是滚子x轴方向的自转角速度
% w2 是内圈 的角速度
n=30;Dw=0.015;Dm=(0.26+0.34)/2;lenroller=0.015;lenrollerline=0.015;arcroller=0.002;CR=0.0001;
deltar0=15e-5; Rp=0.013; %兜孔半径
a0=datafromvb(9)*pi/180;
ballden=7850;e1=2.145e11;e2=2.145e11;e3=2.145e11;o1=0.2808;o2=0.2808;o3=0.2808; 
W1=0*pi/30;W2=9900*pi/30; Fyy=0; Fzz=25000+0.1; Myy=0+0.1; Mzz=0+0.1; 
cucao1=1e-07;cucao2=1e-07;cucao3=1e-07; 
%润滑油基本参数
oilden=885; sita0=20;K=0.14; %常温下的导热系数
niandu0=0.0318; nianya0=1.25e-08; beita0=0.032;   %粘温系数（近似认为不变，wys论文中）   所有参数均为常温下的参数！！！
yindao=1; yindaojianxi=6e-04;   % 引导方式：引导间隙
Dr1=Dm+Dw+deltar0;Dr2=Dm-Dw-deltar0;  % 滚道直径
gama=Dw/Dm;
leixing=5
D=Dw*1000; Dm=Dm*1000; l=(lenroller-2*CR)*1000; aa=0; w2=W2;

if (leixing==1)%深沟球面滚子轴承
bm=1.15  ;
elseif(leixing==2)%向心圆柱滚子轴承
bm=1.1;
elseif(leixing==3)%向心圆锥滚子轴承
bm=1.1 ;
elseif(leixing==4)%实心套圈向心滚针轴承
bm=1.1 
elseif(leixing==5)%冲压滚针轴承
bm=1
elseif(leixing==6)%推力圆锥滚子轴承
bm=1.1 
elseif(leixing==7)%推力调心滚子轴承
bm=1.15
elseif(leixing==8)%推力圆柱滚子轴承
bm=1
elseif(leixing==9)%推力滚针轴承
bm=1 
end


e=9/8;
aa=0;
a=pi.*aa/180;
R=D*cos(a)/Dm
Z=size(ROLLER1,2)%滚子数
m=size(ROLLER1,1)%分段数
%内圈
Qgc1=bm*552*((1-R)^(29/27)/(1+R)^(1/4))*(R/cos(a))^(2/9)*D^(29/27)*(l)^(7/9)*Z^(-1/4)/150
Lu=((Qgc1)./ROLLER1).^4
[m,n]=size(Lu);
G=m*n
%disp('上册方法')
%Lu=Lu.^(-4/4)
%L1=((sum(Lu(:)))/G)^(-4/4)
disp('下册方法')
Lu=Lu.^(-e);
L11=((sum(Lu(:)))/G)^(-1/e);

%外圈
Qgc2=bm*552*((1+R)^(29/27)/(1-R)^(1/4))*(R/cos(a))^(2/9)*D^(29/27)*(l)^(7/9)*Z^(-1/4)/150;
Lv=((Qgc2)./ROLLER2).^(4.5);
Lv=Lv.^(-e);
L2=((sum(Lv(:)))/G)^(-1/e);
Lxia=((sum(Lu(:)))/G+(sum(Lv(:)))/G)^(-1/e)%滚道寿命
%disp('上册方法')
%L=(L11^(-e)+L2^(-e))^(-1/e)
%L=((sum(Lu(:))+sum(Lv(:)))/14)^(-1/e)
%————————————————————————滚子寿命———————————————————————

%内圈滚子
Qgnjk1=464*(1-R)^1.324*(l/m)^(7/9)*(D^(29/27))/(cos(a))^(2/9);
Lgu=((Qgnjk1)./ROLLER1).^(4);
n=w1/w2;
n=sum(n)/length(n);
Lgu=Lgu./n ; %转换成内圈转数
Lgu=Lgu.^(-e)  ;
Lgnei=((sum(Lgu(:)))/G)^(-1/e)
%外圈滚子
Qgnjk2=464*(1+R)^1.324*(l/m)^(7/9)*D^(29/27)/(cos(a))^(2/9);
Lgv=((Qgnjk2)./ROLLER2).^(4.5);
Lgv=Lgv/n;   %转换成内圈转数
Lgv=Lgv.^(-e);
Lgwai=((sum(Lgv(:)))/G)^(-1/e)
L=((sum(Lu(:)))/G+(sum(Lv(:)))/G+(sum(Lgu(:)))/G+(sum(Lgv(:)))/G)^(-1/e)
save life_roller Lxia L
%----------------------------------------打开文件继续写入--------------------------
