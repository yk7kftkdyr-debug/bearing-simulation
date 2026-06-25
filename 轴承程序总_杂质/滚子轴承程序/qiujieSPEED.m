function qiujieSPEED(datafromvb )  % speed11

load Q1;load Q2;load loadi;
% 数据输入部分！！
n=datafromvb(1);Dw=datafromvb(2);Dm=datafromvb(3);lenroller=datafromvb(4);lenrollerline=datafromvb(5);arcroller=datafromvb(6);
deltar0=datafromvb(7); Rp=datafromvb(8); datafromvb(75)=14;%兜孔半径
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

markloadi=1;              % 将滚子按 承载连续性排序，如 12 13 14 1 2 3,  而不是 1 2 3 12 13 14
if loadj>1
for i=1:loadj
Wo(i)=(W2*(1-gama)-W1*(1+gama))/2; 
Wx(i)=(W2-Wo(i))*(1-gama)/gama;
if i==1
    if (loadi(2)-loadi(1))~=1
       markloadi=2;
    end
else
   if (loadi(i)-loadi(i-1))~=1
       markloadi=i;
   end
end
end
end
loadii=[loadi(markloadi:end)  loadi(1:markloadi-1)];  save loadii;    % 将滚子按 承载连续性排序
Ph1ii=[Ph1(markloadi:end)  Ph1(1:markloadi-1)];  save Ph1ii;
Ph2ii=[Ph2(markloadi:end)  Ph2(1:markloadi-1)];  save Ph2ii;
Q1ii=[Q1(markloadi:end)  Q1(1:markloadi-1)];  save Q1ii;
Q2ii=[Q2(markloadi:end)  Q2(1:markloadi-1)];  save Q2ii;

Wononload=(W2*(1-gama)-W1*(1+gama))/2;      Wxnonload=(W2-Wononload)*(1-gama)/gama;
SPEED1=[ Wo Wx Wononload Wxnonload]';
www2=SPEED1; 

for i=1:loadj
    if Q1ii(i)==max(Q1ii)
       markQ1=i;save markQ1
    end
end

nn2=1; speedchange=1e-5; result222=1e4;min2=1e50; wwmin2=www2;
for i=1:(2*loadj+1)
    ddd(i)=0;
end
while( result222>0.001 )
    
    ffSPEED(datafromvb,www2);     load result222;   aaa=result222;
    
    if result222<min2
    wwmin2=www2; save wwmin2;min2=result222
    end
    
    q2=www2;    
    for ij=1:loadj
    q2(ij)=www2(ij)*(1+speedchange);  ffSPEED(datafromvb,q2);load result222;bbb=result222;   ddd(ij)=(bbb-aaa)/(speedchange*www2(ij)); q2=www2;
    q2(loadj+ij)=www2(loadj+ij)*(1+speedchange);   ffSPEED(datafromvb,q2);load result222;bbb=result222;    ddd(loadj+ij)=(bbb-aaa)/(speedchange*www2(loadj+ij));  q2=www2;
    end
    q2(2*loadj+1)=www2(2*loadj+1)*(1+speedchange);   ffSPEED(datafromvb,q2);load result222;bbb=result222;    ddd(2*loadj+1)=(bbb-aaa)/(speedchange*www2(2*loadj+1));  q2=www2;
    q2(2*loadj+2)=www2(2*loadj+2)*(1+speedchange);   ffSPEED(datafromvb,q2);load result222;bbb=result222;    ddd(2*loadj+2)=(bbb-aaa)/(speedchange*www2(2*loadj+2));  q2=www2;
      
    dd=0;
    for j=1:(2*loadj+2)
        dd=dd+ddd(j)^2;
    end

    suanzi=aaa/dd;
  
    qqq=www2
    for k=1:2*loadj+2
    qqq(k)=www2(k)-(suanzi*ddd(k))/1; 
    end
    www2=qqq'
    save www2;      nn2=nn2+1
    if nn2>10
        break
    end
end

save ddd; save wwmin2;



