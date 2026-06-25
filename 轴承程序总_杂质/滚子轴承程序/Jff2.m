function Jff2(datafromvb,loadi, ttt)       %求 滚子和内圈 载荷平衡方程组的 雅可比矩阵！

pp=ttt;  
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
longaxis=Dr1/2+tuoyuandu;  shortaxis=Dr1/2-tuoyuandu;

m=ballden*Dw*Dw*lenroller*pi/4;  %滚子的质量计算
J=m*Dw*Dw/8;             %转动惯量计算
E1=2/(((1-o1^2)/e1)+((1-o3^2)/e3));             %e1，e3为外圈和滚子的弹性模量，o1，o3为外圈，滚子的泊松比
E2=2/(((1-o2^2)/e2)+((1-o3^2)/e3));   %e1，e3为内圈和滚子的弹性模量，o1，o3为内圈，滚子的泊松比
Rr1=Dr1/2*Dw/2/(Dr1/2-Dw/2);  % 滚子与外圈的当量曲率半径！！！
Rr2=Dr2/2*Dw/2/(Dr2/2+Dw/2);  % 滚子与内圈的当量曲率半径！！！

thickness=lenroller/150;
UU1=0.5*0.5*Dm*((1+gama)*(W1+(W2*(1-gama)-W1*(1+gama))/2)+gama*(W2-(W2*(1-gama)-W1*(1+gama))/2)*(1-gama)/gama);
UU2=0.5*0.5*Dm*((1-gama)*(W2-(W2*(1-gama)-W1*(1+gama))/2)+gama*(W2-(W2*(1-gama)-W1*(1+gama))/2)*(1-gama)/gama);


%初始的假设值！！！

z6=0;z7=0;z8=0;z9=0;
loadj=length(loadi);
for i=1:loadj
Wo(i)=(W2*(1-gama)-W1*(1+gama))/2 ;
Wx(i)=(W2-Wo(i))*(1-gama)/gama;
end
Z2=pp(2*loadj+1);Y2=pp(2*loadj+2);sitay=pp(2*loadj+3);sitaz=pp(2*loadj+4);
for i=1:loadj
    z2(i)=0;z4(i)=0;
    afa(i)=0;
    beita(i)=pp(loadj+i);
    delta1(i)=pp(i);
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
for i=1:loadj              %包括了0号滚子，
     sita(i)=2*pi*(loadi(i))/n;   
     Q1e(i)=0;Q2e(i)=0; FP1e(i)=0;FP2e(i)=0; eQ2fenzi(i)=0; eQ2fenmu(i)=0;  
     numcontact11(i)=0;  numcontact22(i)=0;     %   统计承载的片数！！
end     
  for i=1:loadj
    for j=1:loadj
     z2afa(i,j)=0; z2beita(i,j)=0;   z2delta1(i,j)=0;     z2Y2(i)=0;  z2Z2(i)=0;  z2sitay(i)=0;  z2sitaz(i)=0;  
     z4afa(i,j)=0; z4beita(i,j)=0;   z4delta1(i,j)=0;     z4Y2(i)=0;  z4Z2(i)=0;  z4sitay(i)=0;  z4sitaz(i)=0; 
    end      
end
for j=1:n
  z6afa(i)=0; z6beita(j)=0; z6delta1(j)=0; 
  z7afa(i)=0; z7beita(j)=0; z7delta1(j)=0;
  z8afa(i)=0; z8beita(j)=0; z8delta1(j)=0; 
  z9afa(i)=0; z9beita(j)=0; z9delta1(j)=0;
end
z6Y2=0; z6Z2=0; z6sitay=0; z6sitaz=0; 
z7Y2=0; z7Z2=0; z7sitay=0; z7sitaz=0; 
z8Y2=0; z8Z2=0; z8sitay=0; z8sitaz=0; 
z9Y2=0; z9Z2=0; z9sitay=0; z9sitaz=0; 



for i=1:loadj             %包括了0号滚子，
         sita(i)=2*pi*(loadi(i))/n;   
         Q1(i)=0;   Q2(i)=0;     Q1e(i)=0;   Q2e(i)=0;     eQ2fenzi(i)=0;          eQ2fenmu(i)=0;
         Q1afa(i,i)= 0;          Q1beita(i,i)= 0;          Q1delta1(i,i)=0;
         Q1Y2(i)=0;              Q1Z2(i)=0;                Q1sitay(i)=0;           Q1sitaz(i)= 0;
         Q2afa(i,i)=0;           Q2beita(i,i)= 0;          Q2delta1(i,i)= 0;
         Q2Y2(i)=0;              Q2Z2(i)=0;                Q2sitay(i)= 0;          Q2sitaz(i)= 0;  
         Q1eafa(i,i)= 0;          Q1ebeita(i,i)= 0;          Q1edelta1(i,i)=0;
         Q1eY2(i)=0;              Q1eZ2(i)=0;                Q1esitay(i)=0;           Q1esitaz(i)= 0;
         Q2eafa(i,i)=0;           Q2ebeita(i,i)= 0;          Q2edelta1(i,i)= 0;
         Q2eY2(i)=0;              Q2eZ2(i)=0;                Q2esitay(i)= 0;          Q2esitaz(i)= 0;       
     
        eQ2fenziafa(i,i)=0;        eQ2fenzibeita(i,i)=0;        eQ2fenzidelta1(i,i)=0;
        eQ2fenziY2(i)=0;           eQ2fenziZ2(i)=0;        eQ2fenzisitay(i)=0;        eQ2fenzisitaz(i)=0;
        eQ2fenmuafa(i,i)=0 ;        eQ2fenmubeita(i,i)=0;        eQ2fenmudelta1(i,i)=0;
        eQ2fenmuY2(i)=0;           eQ2fenmuZ2(i)=0;        eQ2fenmusitay(i)=0 ;        eQ2fenmusitaz(i)=0;
         
         eQ2afa(i,i)=0;           eQ2beita(i,i)= 0;          eQ2delta1(i,i)= 0;
         eQ2Y2(i)=0;              eQ2Z2(i)=0;                eQ2sitay(i)= 0;          eQ2sitaz(i)= 0;          
       
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
%*****************************************************************************************   
 
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
        %delta2k(i)=delta2(i)-ck+delta2kafa+delta2kbeita-oilh2(i)+debris_shift ;  
%*****************************************************************************************  

%********************************未考虑杂质影响的接触变形*********************************************************  
         delta1k(i)=delta1(i)-ck+delta1kafa+delta1kbeita-oilh1(i)+debris_shift ;    
         delta2k(i)=delta2(i)-ck+delta2kafa+delta2kbeita-oilh2(i)+debris_shift ; 
%*****************************************************************************************  
                
         if delta1k(i)<0
            Q1k=0;
         else
             Q1k= (delta1k(i)*(lenroller*1000)^0.8/3.84e-8)^1.111;             
         end     
        
         if delta2k(i)<0
            Q2k=0;
         else
             Q2k= (delta2k(i)*(lenroller*1000)^0.8/3.84e-8)^1.111;         
         end  
         
         Q1(i)=Q1(i)+Q1k/150;
         Q2(i)=Q2(i)+Q2k/150;
        
        eQ2fenzi(i)=eQ2fenzi(i)+Q2k/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenmu(i)=eQ2fenmu(i)+Q2k/lenroller;
        
            Q1e(i)=Q1e(i)+Q1k*(lenroller/2-(k-0.5)*thickness);
            Q2e(i)=Q2e(i)+Q2k*(lenroller/2-(k-0.5)*thickness);
    
        
      %求解偏导数！！
        if cos(sita(i))>0
             delta1kbeitabeita(i,i)=( -lenroller/2+(k-0.5)*lenroller/150  + Dr2/2*tan(sitay*cos(sita(i))/2) )/cos(beita(i))^2;
             delta2kbeitabeita(i,i)=-( -lenroller/2+(k-0.5)*lenroller/150  + Dr2/2*tan(sitay*cos(sita(i))/2) )/cos(sitay*cos(sita(i))-beita(i))^2;
             delta1kbeitasitay(i)=Dr2/2*cos(sita(i))/2/cos(sitay*cos(sita(i))/2)^2*tan(beita(i));
             delta2kbeitasitay(i)=Dr2/2*cos(sita(i))/2/cos(sitay*cos(sita(i))/2)^2*tan(sitay*cos(sita(i))-beita(i))+( -lenroller/2+(k-0.5)*lenroller/150 + Dr2/2*tan(sitay*cos(sita(i))/2) )*cos(sita(i))/cos(sitay*cos(sita(i))-beita(i))^2;
        else
             delta1kbeitabeita(i,i)=-( -lenroller/2+(k-0.5)*lenroller/150  + Dr2/2*tan(sitay*cos(sita(i))/2) )/cos(beita(i))^2;
             delta2kbeitabeita(i,i)=( -lenroller/2+(k-0.5)*lenroller/150  + Dr2/2*tan(sitay*cos(sita(i))/2) )/cos(sitay*cos(sita(i))+beita(i))^2;
             delta1kbeitasitay(i)=Dr2/2*cos(sita(i))/2/cos(sitay*cos(sita(i))/2)^2*tan(-beita(i));
             delta2kbeitasitay(i)=Dr2/2*cos(sita(i))/2/cos(sitay*cos(sita(i))/2)^2*tan(sitay*cos(sita(i))+beita(i))+( -lenroller/2+(k-0.5)*lenroller/150 + Dr2/2*tan(sitay*cos(sita(i))/2) )*cos(sita(i))/cos(sitay*cos(sita(i))+beita(i))^2;
        end
         delta1kafaafa(i,i)=((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i)) )*(-lenroller/2+(k-0.5)*lenroller/150  )*cos(sitaz*sin(sita(i))-afa(i))/sqrt((Dr1/2)^2-((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i)) )^2   ) ;
         delta2kafaafa(i,i)=((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i)) )*(-lenroller/2+(k-0.5)*lenroller/150  )*cos(sitaz*sin(sita(i))-afa(i))/sqrt((Dr2/2)^2-((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i)) )^2   ) ;
         delta1kafasitaz(i)=-sin(sita(i))*((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i))*(-lenroller/2+(k-0.5)*lenroller/150  )*cos(sitaz*sin(sita(i))-afa(i)) )/sqrt((Dr1/2)^2-((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i)) )^2   ) ;
         delta2kafasitaz(i)=-sin(sita(i))*((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i))*(+-lenroller/2+(k-0.5)*lenroller/150 )*cos(sitaz*sin(sita(i))-afa(i)) )/sqrt((Dr2/2)^2-((-lenroller/2+(k-0.5)*lenroller/150  )*sin(sitaz*sin(sita(i))-afa(i)) )^2   ) ;

         delta1kdelta1(i,i)=1;
         delta1kafa(i,i)=delta1kafaafa(i,i);
         delta1kbeita(i,i)= delta1kbeitabeita(i,i);     
         delta1ksitay(i)=  delta1kbeitasitay(i);        
         delta1ksitaz(i)=delta1kafasitaz(i);
         delta2kdelta1(i,i)=-1;
         delta2kafa(i,i)=delta2kafaafa(i,i);
         delta2kbeita(i,i)=   delta2kbeitabeita(i,i);   
         delta2kY2(i)=sin(sita(i));
         delta2kZ2(i)=cos(sita(i));
         delta2ksitay(i)=   delta2kbeitasitay(i); 
         delta2ksitaz(i)=delta2kafasitaz(i);

         if delta1k(i)<0
                 Q1kafa =0;
                 Q1kbeita =0;
                 Q1kdelta1 =0;
                 Q1kY2 =0;   Q1kZ2 =0;
                 Q1ksitay =0;
                 Q1ksitaz =0;
         else
                 Q1kafa =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta1k(i)^0.111*delta1kafa(i,i);
                 Q1kbeita = ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta1k(i)^0.111*delta1kbeita(i,i);
                 Q1kdelta1 =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta1k(i)^0.111*delta1kdelta1(i,i);
                 Q1kY2 =0;   Q1kZ2 =0;
                 Q1ksitay =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta1k(i)^0.111*delta1ksitay(i);
                 Q1ksitaz =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta1k(i)^0.111*delta1ksitaz(i);
         end
         
         
         if delta2k(i)<0
                 Q2kafa =0;
                 Q2kbeita =0;
                 Q2kdelta1 =0;
                 Q2kY2 =0;   Q2kZ2 =0;
                 Q2ksitay =0;
                 Q2ksitaz =0;
         else
         Q2kafa =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2kafa(i,i);
         Q2kbeita =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2kbeita(i,i);
         Q2kdelta1 = ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2kdelta1(i,i);
         Q2kY2 =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2kY2(i);
         Q2kZ2 =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2kZ2(i);
         Q2ksitay =  ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2ksitay(i);
         Q2ksitaz = ((lenroller*1000)^0.8/3.84e-8)^1.111*1.111*delta2k(i)^0.111*delta2ksitaz(i);         
         end

         
         Q1afa(i,i)= Q1afa(i,i)+Q1kafa/150 ;
         Q1beita(i,i)= Q1beita(i,i)+Q1kbeita/150 ;
         Q1delta1(i,i)= Q1delta1(i,i)+Q1kdelta1/150 ;
         Q1Y2(i)=0; 
         Q1Z2(i)=0;
         Q1sitay(i)= Q1sitay(i)+Q1ksitay/150 ;
         Q1sitaz(i)= Q1sitaz(i)+Q1ksitaz/150 ;
         Q2afa(i,i)= Q2afa(i,i)+Q2kafa/150 ;
         Q2beita(i,i)= Q2beita(i,i)+Q2kbeita/150 ;
         Q2delta1(i,i)= Q2delta1(i,i)+Q2kdelta1/150 ;
         Q2Y2(i)=Q2Y2(i)+Q2kY2/150 ; 
         Q2Z2(i)=Q2Z2(i)+Q2kZ2/150 ;
         Q2sitay(i)= Q2sitay(i)+Q2ksitay/150 ;
         Q2sitaz(i)= Q2sitaz(i)+Q2ksitaz/150 ;   

         
        eQ2fenziafa(i,i)=eQ2fenziafa(i,i)+Q2kafa/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenzibeita(i,i)=eQ2fenzibeita(i,i)+Q2kbeita/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenzidelta1(i,i)=eQ2fenzidelta1(i,i)+Q2kdelta1/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenziY2(i)=eQ2fenziY2(i)+Q2kY2/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenziZ2(i)=eQ2fenziZ2(i)+Q2kZ2/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenzisitay(i)=eQ2fenzisitay(i)+Q2ksitay/lenroller*(lenroller/2-(k-0.5)*thickness);
        eQ2fenzisitaz(i)=eQ2fenzisitaz(i)+Q2ksitaz/lenroller*(lenroller/2-(k-0.5)*thickness);
        
        eQ2fenmuafa(i,i)=eQ2fenmuafa(i,i)+Q2kafa/lenroller ;
        eQ2fenmubeita(i,i)=eQ2fenmubeita(i,i)+Q2kbeita/lenroller ;
        eQ2fenmudelta1(i,i)=eQ2fenmudelta1(i,i)+Q2kdelta1/lenroller ;
        eQ2fenmuY2(i)=eQ2fenmuY2(i)+Q2kY2/lenroller ;
        eQ2fenmuZ2(i)=eQ2fenmuZ2(i)+Q2kZ2/lenroller ;
        eQ2fenmusitay(i)=eQ2fenmusitay(i)+Q2ksitay/lenroller ;
        eQ2fenmusitaz(i)=eQ2fenmusitaz(i)+Q2ksitaz/lenroller ;
        

            Q1eafa(i,i)=Q1eafa(i,i)+Q1kafa*(lenroller/2-(k-0.5)*thickness);
            Q1ebeita(i,i)=Q1ebeita(i,i)+Q1kbeita*(lenroller/2-(k-0.5)*thickness);
            Q1edelta1(i,i)=Q1edelta1(i,i)+Q1kdelta1*(lenroller/2-(k-0.5)*thickness);
            Q1eY2(i)=Q1eY2(i)+Q1kY2*(lenroller/2-(k-0.5)*thickness);
            Q1eZ2(i)=Q1eZ2(i)+Q1kZ2*(lenroller/2-(k-0.5)*thickness);
            Q1esitay(i)=Q1esitay(i)+Q1ksitay*(lenroller/2-(k-0.5)*thickness);
            Q1esitaz(i)=Q1esitaz(i)+Q1ksitaz*(lenroller/2-(k-0.5)*thickness);
            
            Q2eafa(i,i)=Q2eafa(i,i)+Q2kafa*(lenroller/2-(k-0.5)*thickness);
            Q2ebeita(i,i)=Q2ebeita(i,i)+Q2kbeita*(lenroller/2-(k-0.5)*thickness);
            Q2edelta1(i,i)=Q2edelta1(i,i)+Q2kdelta1*(lenroller/2-(k-0.5)*thickness);
            Q2eY2(i)=Q2eY2(i)+Q2kY2*(lenroller/2-(k-0.5)*thickness);
            Q2eZ2(i)=Q2eZ2(i)+Q2kZ2*(lenroller/2-(k-0.5)*thickness);
            Q2esitay(i)=Q2esitay(i)+Q2ksitay*(lenroller/2-(k-0.5)*thickness);
            Q2esitaz(i)=Q2esitaz(i)+Q2ksitaz*(lenroller/2-(k-0.5)*thickness);            
        
    end     % k循环结束！！
    
    
    if eQ2fenmu(i)<=0                               %代表第i个滚子没有承载 
        eQ2(i)=0;                                   % 为保证计算中间过程不出错，将eQ2置0
        eQ2afa(i,i)=0;        eQ2beita(i,i)=0;
        eQ2delta1(i,i)=0;        eQ2Y2(i)=0;
        eQ2Z2(i)=0;        eQ2sitay(i)=0;        eQ2sitaz(i)=0;
       
    else
        eQ2(i)=eQ2fenzi(i)/eQ2fenmu(i);
    
       eQ2afa(i,i)=(eQ2fenziafa(i,i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmuafa(i))/eQ2fenmu(i)^2;
       eQ2beita(i,i)=(eQ2fenzibeita(i,i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmubeita(i))/eQ2fenmu(i)^2;
       eQ2delta1(i,i)=(eQ2fenzidelta1(i,i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmudelta1(i))/eQ2fenmu(i)^2;
       eQ2Y2(i)=(eQ2fenziY2(i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmuY2(i))/eQ2fenmu(i)^2;
       eQ2Z2(i)=(eQ2fenziZ2(i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmuZ2(i))/eQ2fenmu(i)^2;
       eQ2sitay(i)=(eQ2fenzisitay(i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmusitay(i))/eQ2fenmu(i)^2;
       eQ2sitaz(i)=(eQ2fenzisitaz(i)*eQ2fenmu(i)-eQ2fenzi(i)*eQ2fenmusitaz(i))/eQ2fenmu(i)^2;
    end
        
    Fc(i)=0.5*m*Dm*Wo(i)^2;
    
    z2(i)=Q2(i)-Q1(i)+Fc(i);
    z4(i)=Q2e(i)-Q1e(i);    
        z6=z6+Q2(i)*cos(sita(i));
        z7=z7+Q2(i)*sin(sita(i));
        z8=z8+Q2(i)*eQ2(i)*cos(sita(i));
        z9=z9+Q2(i)*eQ2(i)*sin(sita(i));
  


  %求非线性方程组的偏导数！

      z2afa(i,i)=Q2afa(i,i)-Q1afa(i,i);
      z2beita(i,i)=Q2beita(i,i)-Q1beita(i,i);
      z2delta1(i,i)=Q2delta1(i,i)-Q1delta1(i,i);
      z2Y2(i)=Q2Y2(i)-Q1Y2(i);
      z2Z2(i)=Q2Z2(i)-Q1Z2(i);
      z2sitay(i)=Q2sitay(i)-Q1sitay(i);
      z2sitaz(i)=Q2sitaz(i)-Q1sitaz(i);
      
      z4afa(i,i)=Q2eafa(i,i)-Q1eafa(i,i);
      z4beita(i,i)=Q2ebeita(i,i)-Q1ebeita(i,i);
      z4delta1(i,i)=Q2edelta1(i,i)-Q1edelta1(i,i);
      z4Y2(i)=Q2eY2(i)-Q1eY2(i);
      z4Z2(i)=Q2eZ2(i)-Q1eZ2(i);
      z4sitay(i)=Q2esitay(i)-Q1esitay(i);
      z4sitaz(i)=Q2esitaz(i)-Q1esitaz(i);
      
        z6afa(i)=z6afa(i)+Q2afa(i,i)*cos(sita(i));
        z6beita(i)=z6beita(i)+Q2beita(i,i)*cos(sita(i));
        z6delta1(i)=z6delta1(i)+Q2delta1(i,i)*cos(sita(i));
        z6Y2=z6Y2+Q2Y2(i)*cos(sita(i));
        z6Z2=z6Z2+Q2Z2(i)*cos(sita(i));
        z6sitay=z6sitay+Q2sitay(i)*cos(sita(i));
        z6sitaz=z6sitaz+Q2sitaz(i)*cos(sita(i));
      
        z7afa(i)=z7afa(i)+Q2afa(i)*sin(sita(i));
        z7beita(i)=z7beita(i)+Q2beita(i,i)*sin(sita(i));
        z7delta1(i)=z7delta1(i)+Q2delta1(i,i)*sin(sita(i));
        z7Y2=z7Y2+Q2Y2(i)*sin(sita(i));
        z7Z2=z7Z2+Q2Z2(i)*sin(sita(i));
        z7sitay=z7sitay+Q2sitay(i)*sin(sita(i));
        z7sitaz=z7sitaz+Q2sitaz(i)*sin(sita(i));
       
        z8afa(i)=z8afa(i)+Q2afa(i,i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2afa(i,i)*cos(sita(i));
        z8beita(i)=z8beita(i)+Q2beita(i,i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2beita(i,i)*cos(sita(i));
        z8delta1(i)=z8delta1(i)+Q2delta1(i,i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2delta1(i,i)*cos(sita(i));
        z8Y2=z8Y2+Q2Y2(i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2Y2(i)*cos(sita(i));
        z8Z2=z8Z2+Q2Z2(i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2Z2(i)*cos(sita(i));
        z8sitay=z8sitay+Q2sitay(i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2sitay(i)*cos(sita(i));
        z8sitaz=z8sitaz+Q2sitaz(i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2sitaz(i)*cos(sita(i));
        
        z9afa(i)=z9afa(i)+Q2afa(i,i)*eQ2(i)*cos(sita(i))+Q2(i)*eQ2afa(i,i)*sin(sita(i));
        z9beita(i)=z9beita(i)+Q2beita(i,i)*eQ2(i)*sin(sita(i))+Q2(i)*eQ2beita(i,i)*sin(sita(i));
        z9delta1(i)=z9delta1(i)+Q2delta1(i,i)*eQ2(i)*sin(sita(i))+Q2(i)*eQ2delta1(i,i)*sin(sita(i));
        z9Y2=z9Y2+Q2Y2(i)*eQ2(i)*sin(sita(i))+Q2(i)*eQ2Y2(i)*sin(sita(i));
        z9Z2=z9Z2+Q2Z2(i)*eQ2(i)*sin(sita(i))+Q2(i)*eQ2Z2(i)*sin(sita(i));
        z9sitay=z9sitay+Q2sitay(i)*eQ2(i)*sin(sita(i))+Q2(i)*eQ2sitay(i)*sin(sita(i));
        z9sitaz=z9sitaz+Q2sitaz(i)*eQ2(i)*sin(sita(i))+Q2(i)*eQ2sitaz(i)*sin(sita(i));        

          
end   %对应于for的最后一个end


% 构建雅可比矩阵！！！！
Jzz=zeros(2*loadj+3,2*loadj+3);

for i=1:loadj
    for j=1:loadj

    Jzz(i,j)=z2delta1(i,j);                         Jzz(i,(loadj+j))=z2beita(i,j);                Jzz(i,(2*loadj+1))=z2Z2(i);  
    Jzz(i,(2*loadj+2))=z2Y2(i);                     Jzz(i,(2*loadj+3))=z2sitay(i);                
   
    Jzz(loadj+i,j)=z4delta1(i,j);               Jzz(loadj+i,(loadj+j))=z4beita(i,j);      Jzz(loadj+i,(2*loadj+1))=z4Z2(i);  
    Jzz(loadj+i,(2*loadj+2))=z4Y2(i);           Jzz(loadj+i,(2*loadj+3))=z4sitay(i);      

    end
end

for j=1:loadj
    
   Jzz(2*loadj+1,j)=-z6delta1(j);                Jzz(2*loadj+1,(loadj+j))=-z6beita(j);       Jzz(2*loadj+1,(2*loadj+1))=-z6Z2  ;     
   Jzz(2*loadj+1,(2*loadj+2))=-z6Y2;              Jzz(2*loadj+1,(2*loadj+3))=-z6sitay;  
   
    
    Jzz(2*loadj+2,j)=-z7delta1(j);               Jzz(2*loadj+2,(loadj+j))=-z7beita(j);        Jzz(2*loadj+2,(2*loadj+1))=-z7Z2  ;     
    Jzz(2*loadj+2,(2*loadj+2))=-z7Y2;             Jzz(2*loadj+2,(2*loadj+3))=-z7sitay;   
   
   
    Jzz(2*loadj+3,j)= z8delta1(j);               Jzz(2*loadj+3,(loadj+j))= z8beita(j);        Jzz(2*loadj+3,(2*loadj+1))= z8Z2;     
    Jzz(2*loadj+3,(2*loadj+2))=z8Y2;             Jzz(2*loadj+3,(2*loadj+3))= z8sitay;      
  
end

save Jzz






