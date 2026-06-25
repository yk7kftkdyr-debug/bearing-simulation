function ff1(datafromvb )   % 求解 承载滚动体数目！！

% 数据输入部分！！
n=datafromvb(1);Dw=datafromvb(2);Dm=datafromvb(3);lenroller=datafromvb(4);lenrollerline=datafromvb(5);arcroller=datafromvb(6);
deltar0=datafromvb(7); Rp=datafromvb(8);datafromvb(75)=14; %兜孔半径
a0=datafromvb(9)*pi/180;
ballden=datafromvb(59);e1=datafromvb(11);e2=datafromvb(12);e3=datafromvb(13);o1=datafromvb(14);o2=datafromvb(15);o3=datafromvb(16); 
W1=datafromvb(17)*pi/30;W2=datafromvb(18)*pi/30; Fyy=datafromvb(19); Fzz=datafromvb(20)+0.1; Myy=datafromvb(21)+0.1; Mzz=datafromvb(22)+0.1; 
cucao1=datafromvb(23);cucao2=datafromvb(24);cucao3=datafromvb(25); 
%润滑油基本参数
oilden=datafromvb(26); sita0=datafromvb(27);K=datafromvb(28); %常温下的导热系数
niandu0=datafromvb(29); nianya0=datafromvb(30); beita0=datafromvb(31);   %粘温系数（近似认为不变，wys论文中）   所有参数均为常温下的参数！！！
yindao=datafromvb(32); yindaojianxi=datafromvb(33);   % 引导方式：引导间隙
Dr1=Dm+Dw+deltar0;Dr2=Dm-Dw-deltar0;  % 滚道直径
gama=Dw/Dm;

%初始的假设值！！！
for i=1:n
Wo(i)=(W2*(1-gama)-W1*(1+gama))/2 ;
Wx(i)=(W2-Wo(i))*(1-gama)/gama;
end

m=ballden*Dw*Dw*lenroller*pi/4;  %滚子的质量计算
J=m*Dw*Dw/8;             %转动惯量计算
E1=2/(((1-o1^2)/e1)+((1-o3^2)/e3));             %e1，e3为外圈和滚子的弹性模量，o1，o3为外圈，滚子的泊松比
E2=2/(((1-o2^2)/e2)+((1-o3^2)/e3));   %e1，e3为内圈和滚子的弹性模量，o1，o3为内圈，滚子的泊松比

%求静态载荷分布！！（静态时 滚子与外、内圈赫兹接触载荷相同！）

Q0=4.6*sqrt(Fyy^2+Fzz^2)/n; 
Q00=2*Q0;
while(abs(Q0-Q00)>1e-0)            %收敛条件!!
    Q0=Q00;
deltamax=7.66e-8*Q0^0.9/(lenroller*1000)^0.8;
loadcoeff1=deltamax/(2*deltamax+deltar0);
loadangle1=acos(deltar0/(2*deltamax+deltar0));
loadn=floor(n*loadangle1/6.28);
temp1=0;
for i=1:loadn                    % 这里不包括0号滚子.
    sita(i)=2*pi*i/n;  %求各个滚子的方位角  
    temp1=temp1+(1-(1-cos(sita(i)))/(2*loadcoeff1))^(10/9)*cos(sita(i));
end
Q00=sqrt(Fyy^2+Fzz^2)/( 1+2*temp1 );
end
% 求得承载滚子数为  2*loadn+1  个！
save loadn  

  

  