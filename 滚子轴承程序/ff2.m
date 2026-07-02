function ff2(datafromvb,loadi,ttt1 )   % 求解 承载区  滚动体的1  3 两个方程组成的 方程组

pp=ttt1;  
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
Dr1=Dm+Dw+deltar0;Dr2=Dm-Dw-deltar0;  % 滚道直径
gama=Dw/Dm;
% 对于非圆套圈滚道的计算！！
longaxis=Dr1/2 + tuoyuandu;  shortaxis=Dr1/2 - tuoyuandu;

Rr1=Dr1/2*Dw/2/(Dr1/2-Dw/2);  % 滚子与外圈的当量曲率半径！！！
Rr2=Dr2/2*Dw/2/(Dr2/2+Dw/2);  % 滚子与内圈的当量曲率半径！！！
E1=2/(((1-o1^2)/e1)+((1-o3^2)/e3));             %e1，e3为外圈和球的弹性模量，o1，o3为外圈，球的泊松比
E2=2/(((1-o2^2)/e2)+((1-o3^2)/e3));             %e1，e3为内圈和球的弹性模量，o1，o3为内圈，球的泊松比

m=ballden*Dw*Dw*lenroller*pi/4;  %滚子的质量计算
J=m*Dw*Dw/8;             %转动惯量计算
%初始的假设值！！！
z6=0;z7=0;z8=0;z9=0;                      QQ1=[];QQ2=[];
loadj=length(loadi);
for i=1:loadj
Wo(i)=(W2*(1-gama)-W1*(1+gama))/2; 
Wx(i)=(W2-Wo(i))*(1-gama)/gama;
end
Z2=pp(2*loadj+1);Y2=pp(2*loadj+2);sitay=pp(2*loadj+3);sitaz=pp(2*loadj+4);
for i=1:loadj
    z2(i)=0;z4(i)=0;  afa(i)=0;  
    delta1(i)=pp(i);  beita(i)=pp(loadj+i);  
    sita(i)=2*pi*(loadi(i))/n; 
    % 考虑非圆滚道计算
     deltachange=sqrt( 1/(longaxis^(-2)+tan(sita(i))^2/shortaxis^2) * (1+tan(sita(i))^2) )-Dr1/2 ;
     deltar=deltar0+deltachange;
    delta2(i)=Y2*sin(sita(i))+Z2*cos(sita(i))-deltar/2-delta1(i);
end
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
     Q1(i)=0;Q2(i)=0;Q1e(i)=0;Q2e(i)=0; eQ2fenzi(i)=0; eQ2fenmu(i)=0;  
     numcontact11(i)=0;  numcontact22(i)=0;     %   统计承载的片数！！
     
   % 最小油膜厚度计算!!!
   %yang解释：等温、考虑热效应、考虑表面纹理参数，使用对应油膜厚度计算公式；考虑杂质影响，在接触变形delta1k，delta2k上加上特征位移ud  
   
%***************************等温条件下膜厚计算公式******************************************      micro_config = load_micro_interface_config();
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
oilh1(i)=thermal_factor1*texture_factor*2.65*Rr1*(nianya0*E1)^0.54*(niandu0*UU1/(E1*Rr1))^0.7*(Q10(i)/(lenroller*E1*Rr1))^(-0.13);oilh2(i)=thermal_factor2*texture_factor*2.65*Rr2*(nianya0*E2)^0.54*(niandu0*UU2/(E2*Rr2))^0.7*(Q20(i)/(lenroller*E2*Rr2))^(-0.13);
 
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
%oilh1(i)=2.65*Rr1*(nianya0*E1)^0.54*(niandu0*UU1/(E1*Rr1))^0.7*(Q10(i)/(lenroller*E1*Rr1))^(-0.13);
%oilh2(i)=2.65*Rr2*(nianya0*E2)^0.54*(niandu0*UU2/(E2*Rr2))^0.7*(Q20(i)/(lenroller*E2*Rr2))^(-0.13);
%*****************************************************************************************   

    
    for k=1:150
        if (k-0.5)*lenroller/150<=(lenroller-lenrollerline)/2
            ck=sqrt(arcroller^2-lenrollerline^2/4)-sqrt( arcroller^2-(lenroller/2 -(k-0.5)*lenroller/150 )^2 );
        elseif (k-0.5)*lenroller/150>=(lenroller+lenrollerline)/2 && (k-0.5)*lenroller/150<=lenroller
             ck=sqrt(arcroller^2-lenrollerline^2/4)-sqrt( arcroller^2-( (k-0.5)*lenroller/150-lenroller/2 )^2 );
        else
            ck=0;
        end

        if cos(sita(i))>0
            delta1kbeita=(  -lenroller/2+(k-0.5)*lenroller/150 + Dr2/2*tan(sitay*cos(sita(i))/2) )*tan(beita(i));
            delta2kbeita=(  -lenroller/2+(k-0.5)*lenroller/150+ Dr2/2*tan(sitay*cos(sita(i))/2) )*tan(sitay*cos(sita(i))-beita(i));
        else
            delta1kbeita=(  -lenroller/2+(k-0.5)*lenroller/150 + Dr2/2*tan(sitay*cos(sita(i))/2) )*tan(-beita(i));
            delta2kbeita=(  -lenroller/2+(k-0.5)*lenroller/150+ Dr2/2*tan(sitay*cos(sita(i))/2) )*tan(sitay*cos(sita(i))+beita(i));
        end
        
       delta1kafa=sqrt((Dr1/2)^2-((-lenroller/2+(k-0.5)*lenroller/150)*sin(sitaz*sin(sita(i))-afa(i)) )^2   )-Dr1/2;
       delta2kafa=sqrt((Dr2/2)^2-((-lenroller/2+(k-0.5)*lenroller/150)*sin(sitaz*sin(sita(i))-afa(i)) )^2   )-Dr2/2;

%*******************************考虑杂质影响的接触变形********************************************************** 
        %delta1k(i)=delta1(i)-ck+delta1kafa+delta1kbeita-oilh1(i)+debris_shift;    
        %delta2k(i)=delta2(i)-ck+delta2kafa+delta2kbeita-oilh2(i)+debris_shift;  
%*****************************************************************************************  

%********************************未考虑杂质影响的接触变形*********************************************************  
         delta1k(i)=delta1(i)-ck+delta1kafa+delta1kbeita-oilh1(i)+debris_shift ;    
         delta2k(i)=delta2(i)-ck+delta2kafa+delta2kbeita-oilh2(i)+debris_shift ; 
%*****************************************************************************************  
        
         if delta1k(i)<0
            Q1k=0;
         else
             numcontact11(i)=numcontact11(i)+1;
             Q1k= (delta1k(i)*(lenroller*1000)^0.8/3.84e-8)^1.111;   
             widthcontact11(k,i)=sqrt(8*Q1k*Rr1/(pi*lenroller*E1));
             Ph11(k,i)=sqrt(Q1k*E1/(2*pi*lenroller*Rr1));
         end     
        
         if delta2k(i)<0
            Q2k=0;
         else
             numcontact22(i)=numcontact22(i)+1;
             Q2k= (delta2k(i)*(lenroller*1000)^0.8/3.84e-8)^1.111;    
             widthcontact22(k,i)=sqrt(8*Q2k*Rr2/(pi*lenroller*E2));
             Ph22(k,i)=sqrt(Q2k*E2/(2*pi*lenroller*Rr2));
         end  
        
        QQ1(k,i)=Q1k/lenroller;   QQ2(k,i)=Q2k/lenroller;                 % 看 载荷强度 分布情况！！！

        Q1(i)=Q1(i)+Q1k/150;   % 简化后的式子，应该是 lenroller/150*Q1k/lenroller;
        Q2(i)=Q2(i)+Q2k/150;
        
        eQ2fenzi(i)=eQ2fenzi(i)+Q2k/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenmu(i)=eQ2fenmu(i)+Q2k/lenroller;
        
            Q1e(i)=Q1e(i)+Q1k*(lenroller/2-(k-0.5)*thickness);
            Q2e(i)=Q2e(i)+Q2k*(lenroller/2-(k-0.5)*thickness);
    end
    
    lengthcontact1(i)=numcontact11(i)*thickness; lengthcontact2(i)=numcontact22(i)*thickness;   % 接触长度计算！！！
    widthcontact1(i)=max( widthcontact11(:,i) ); widthcontact2(i)=max( widthcontact22(:,i) );  
    Ph1(i)=max( Ph11(:,i) ); Ph2(i)=max( Ph22(:,i) );
    
    if eQ2fenmu(i)<=0          %代表第i个滚子没有承载 
        eQ2(i)=0;            % 为保证计算中间过程不出错，将eQ2置0
    else
        eQ2(i)=eQ2fenzi(i)/eQ2fenmu(i);
    end

    Fc(i)=0.5*m*Dm*Wo(i)^2;

z2(i)=Q2(i)-Q1(i)+Fc(i);
z4(i)=Q2e(i)-Q1e(i);    

if yindao==2   %内圈引导
   %求引导面/保持架之间的法向力和摩擦力！
   Kbp=11/yindaojianxi;
   deltayindao=sqrt(Y2^2+Z2^2);     % 偏离原来位置的距离！
   Fyindao=Kbp*abs(deltayindao);
   fyindao=0.05*Fyindao;    
   zuoyongjiao=atan(Y2/(Z2+1));    
   z6=z6+Q2(i)*cos(sita(i))+Fyindao*cos(zuoyongjiao);
   z7=z7+Q2(i)*sin(sita(i))+Fyindao*sin(zuoyongjiao);
else     %外圈引导 引导面对内圈没有作用力
    z6=z6+Q2(i)*cos(sita(i));
    z7=z7+Q2(i)*sin(sita(i));
    Fyindao=0;    fyindao=0;        zuoyongjiao=0; 
end
save Fyindao;save fyindao;save zuoyongjiao

z8=z8+Q2(i)*eQ2(i)*cos(sita(i));
z9=z9+Q2(i)*eQ2(i)*sin(sita(i));
      
end

z6=Fzz-z6;
z7=Fyy-z7;
z8=Myy+z8;
z9=Mzz+z9;

zz1=[ z2';z4';z6;z7;z8; ];save zz1;
zz2=[ z2';z4';z7;z6;z8;z9 ];save zz2;    %计算4阶刚度使用！！！

result111=0;
for i=1:2*loadj+3
    result111=result111+zz1(i)^2;
end
save result111

save Q1;save Q2;        %  接触载荷  ！！！
save QQ1;save QQ2;  save Ph11;save Ph22    %  载荷强度 分布情况！！！
save lengthcontact1;save widthcontact1;save Ph1;   %接触长度、接触宽度、接触应力（都是计算的各个圆片上的最大值）
save lengthcontact2;save widthcontact2;save Ph2;
                  
              