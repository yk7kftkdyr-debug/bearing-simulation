function qiujieLOAD3(datafromvb)                %牛顿 法求解非线性方程组 得到承载区载荷分布

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

qiujieLOAD(datafromvb);
load www1;   
%初值给定
pp=[];pp(3*loadj+1 )=www1(2*loadj+1); pp(3*loadj+2)=www1(2*loadj+2);pp(3*loadj+3)=www1(2*loadj+3);pp(3*loadj+4)=Mzz/1e2;
for i=1:loadj
    pp(i)=www1(i);
    pp(loadj+i )=www1(loadj+i);  % beita(i)=0
    pp(2*loadj+i )=Mzz/1e3;  % afa(i)=0
end  
% 求解
www12=pp';result112=1e4; nn2=1; temp11=[];  afachange=1e-5; min2=1e50;
for i=1:(3*loadj+4)
    ddd(i)=0;
end
while( result112>0.001 )
    
    ff3(datafromvb,loadi,www12);   load result112;   aaa=result112;
    if result112<min2
    wwmin2=www12; save wwmin2;min2=result112
    end
    q2=www12;    
    for ij=1:loadj
    q2(2*loadj+ij)=www12(2*loadj+ij)*(1+afachange);   ff3(datafromvb,loadi,q2);  load result112;bbb=result112;    ddd(2*loadj+ij)=(bbb-aaa)/(afachange*www12(2*loadj+ij));  q2=www12;
    end
    q2(3*loadj+4)=www12(3*loadj+4)*(1+afachange);     ff3(datafromvb,loadi,q2);  load result112;bbb=result112;    ddd(3*loadj+4)=(bbb-aaa)/(afachange*www12(3*loadj+4));  q2=www12;
    dd=0;
    for j=1:(3*loadj+4)
        dd=dd+ddd(j)^2;
    end
    suanzi=aaa/dd;
    qqq=www12;
    for k=1:3*loadj+4
    qqq(k)=www12(k)-(suanzi*ddd(k)); 
    end
    www12=qqq'
    save www12;      nn2=nn2+1
    if nn2>20
        break
    end
end

www12=wwmin2;  
ff3(datafromvb,loadi,www12); 
save www12



