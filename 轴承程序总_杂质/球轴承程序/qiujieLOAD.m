function qiujieLOAD(datafromvb)   %牛顿 法求解非线性方程组


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

deltac0=Rp/2-Dw/2;
deltar0=2*(1-cos(a0))*(f1+f2-1)*Dw;               %求原始径向间隙！
deltaa0=2*sin(a0)*(f1+f2-1)*Dw ;                 %求原始轴向间隙！
Dy1=Dm+Dw+deltar0;Dy2=Dm-Dw-deltar0;  % 沟底直径
Dr1=Dy1-2*dangbianxishu1*Dw;Dr2=Dy2+2*dangbianxishu2*Dw;  % 挡边直径（引导面直径）  档边高系数给出！！！！

result000=1e4;result111=1e4;result222=1e4;nn=1;  i=1;iii=1;j=1;
min1=1e20;

loadj=n; loadi=1:n; mark1=0;         %mark0 is the mark of computation!!!
% 给定变量初值！！x r X2 Y2  Z2 sitay sitaz初始值
x=[];r=[];
for i=1:loadj
    x=[x Fxx/1e6];   r=[r 1/1e30];
end
X2=Fxx/1e30;  Y2=(Fyy+1)/1e30;Z2=Fzz/1e30;sitay=Myy/1e8;sitaz=Mzz/1e8;       % Z2的初值随着Fxx的变化而变化，且分母也要变化，当Fxx较大时如果分母>1e8,则计算错误！！  
pppp1=[x r X2 Y2  Z2 sitay sitaz]';   %此处x r为相对球坐标系的位移！！！！！！！！！！
www3=pppp1;


ffLOAD(datafromvb,loadi,www3);      load result333; load result444;
nnnum=1;  

if Myy>=1 |  Mzz>=1             % 当存在My Mz时

while (loadj>1 & mark1~=2)         %滚动体承载个数循环
while( result333>0.01 )        %单次计算 循环
    ffLOAD(datafromvb,loadi,www3);      load result333;
    load uu; www3=uu;
    JffLOAD(datafromvb,loadi,www3);    load Jzz;
    www3=[(www3(1:2*loadj+5)-Jzz\zzzz3(1:2*loadj+5))];
    
    for i=1:loadj
        if imag(www3(i))~=0
            mark1=1;break
        end
    end
    if mark1==1
        break
    end
    if result333>1e60
            mark1=1;break
    end
    if mark1==1
        break
    end
    nnnum=nnnum+1;
    if nnnum>60
        mark1=1;    %说明该载荷下 数为 loadj 无法收敛！！！  出错!
    break
    end
end

load q1q2a1a2;
for i=1:loadj
    if   Q1(i)>(sqrt(Fxx^2+Fyy^2+Fzz^2)+Fz(i)) |  Q1(i)<1 | Q2(i)<1 | a1(i)<0 | a2(i)<0       %%判断是否出错！！（即滚动体径向位移是否为负 或者接触角是否<0 ）
        mark1=1;
        break
    end
end
if mark1==1      % 表示计算出错,此时要去掉 承载  滚动体个数!!!!
  if Fyy==0    
    if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[loadi(1:loadj/2-1) loadi(loadj/2+1:loadj)];           % 只去掉一个中间的球
       else
           loadi=[loadi(1:floor(loadj/2))  loadi(floor(loadj/2)+2:loadj)];  % 去掉1个球
    end
  else               %当Fyy~=0 时，不能按照 上方 承载大于下方计算！！！ 应该是 合力 处载荷最大！！ 而不是0处 球！！！
     lipianjiao=atan(Fyy/(Fzz+1));   % Fyy 与 Fzz 的合力与Z 轴的夹角！！！ 
     xishu=lipianjiao/(2*pi);    % 当存在Fyy 时， 去掉承载滚动体个数时 就不能从最下端开始了，而应该沿 合力的反方向！！！
     if (ceil(n*xishu)-(loadj/2-1))<1        % 当可以从1号球开始时
      if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[loadi(1:ceil(n*xishu))  loadi(ceil(n*xishu)+1:ceil(n*xishu)+loadj/2-1)  loadi(ceil(n*xishu)+loadj/2+1:loadj)];   % 只去掉一个中间的球
       else
           loadi=[loadi(1:ceil(n*xishu))  loadi(ceil(n*xishu)+1:ceil(n*xishu)+floor(loadj/2)) loadi(ceil(n*xishu)+floor(loadj/2)+2:loadj)];    % 去掉1个球
       end
     else         % 当 不 可以从1号球开始时
        if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[ (ceil(n*xishu)-(loadj/2-1)):ceil(n*xishu)+loadj/2-1 ];  % 只去掉一个中间的球
       else
           loadi=[  (ceil(n*xishu)-floor(loadj/2)+1):ceil(n*xishu)+floor(loadj/2)   ];    % 去掉1个球
       end 
     end    
  end  
else
    mark1=2;              %说明收敛且没有出错！！！！
    break
end

loadj=length(loadi);   nnnum=1;
% 给定变量初值！！
x=[];r=[];
for i=1:loadj
    x=[x Fxx/1e6];   r=[r 1/1e30];
end
X2=Fxx/1e30;  Y2=(Fyy+1)/1e30;Z2=Fzz/1e30;sitay=Myy/1e8;sitaz=Mzz/1e8;       % Z2的初值随着Fxx的变化而变化，且分母也要变化，当Fxx较大时如果分母>1e8,则计算错误！！  
pppp1=[x r X2 Y2  Z2 sitay sitaz]';   %此处x r为相对球坐标系的位移！！！！！！！！！！
www3=pppp1; 
ffLOAD(datafromvb,loadi,www3);load result333;  ddd1=1e-5;    mark1=0;    % mark1要置0
end

else          % 当 不 存在My  Mz 时
    
while (loadj>1 & mark1~=2)         %滚动体承载个数循环
while( result444>0.01 )        %单次计算 循环
    ffLOAD(datafromvb,loadi,www3);      load result444; 
    load uu; www3=uu;
    JffLOAD(datafromvb,loadi,www3);    load Jzz;
    www3=[(www3(1:2*loadj+3)-Jzz(1:2*loadj+3,1:2*loadj+3)\zzzz3(1:2*loadj+3)) ; Myy/1e8; Mzz/1e8 ];    
    
    for i=1:loadj
        if imag(www3(i))~=0
            mark1=1;break
        end
    end
    if mark1==1
        break
    end
    if result333>1e60
            mark1=1;break
    end
    if mark1==1
        break
    end
    nnnum=nnnum+1;
    if nnnum>60
        mark1=1;    %说明该载荷下 数为 loadj 无法收敛！！！  出错!
    break
    end
end

load q1q2a1a2;
for i=1:loadj
    if   Q1(i)>(sqrt(Fxx^2+Fyy^2+Fzz^2)+Fz(i)) |  Q1(i)<1 | Q2(i)<1 | a1(i)<0 | a2(i)<0       %%判断是否出错！！（即滚动体径向位移是否为负 或者接触角是否<0 ）
        mark1=1;
        break
    end
end
if mark1==1      % 表示计算出错,此时要去掉 承载  滚动体个数!!!!
  if Fyy==0    
    if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[loadi(1:loadj/2-1) loadi(loadj/2+1:loadj)];           % 只去掉一个中间的球
       else
           loadi=[loadi(1:floor(loadj/2))  loadi(floor(loadj/2)+2:loadj)];  % 去掉1个球
    end
  else               %当Fyy~=0 时，不能按照 上方 承载大于下方计算！！！ 应该是 合力 处载荷最大！！ 而不是0处 球！！！
     lipianjiao=atan(Fyy/(Fzz+1));   % Fyy 与 Fzz 的合力与Z 轴的夹角！！！ 
     xishu=lipianjiao/(2*pi);    % 当存在Fyy 时， 去掉承载滚动体个数时 就不能从最下端开始了，而应该沿 合力的反方向！！！
     if (ceil(n*xishu)-(loadj/2-1))<1        % 当可以从1号球开始时
      if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[loadi(1:ceil(n*xishu))  loadi(ceil(n*xishu)+1:ceil(n*xishu)+loadj/2-1)  loadi(ceil(n*xishu)+loadj/2+1:loadj)];   % 只去掉一个中间的球
       else
           loadi=[loadi(1:ceil(n*xishu))  loadi(ceil(n*xishu)+1:ceil(n*xishu)+floor(loadj/2)) loadi(ceil(n*xishu)+floor(loadj/2)+2:loadj)];    % 去掉1个球
       end
     else         % 当 不 可以从1号球开始时
        if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[ (ceil(n*xishu)-(loadj/2-1)):ceil(n*xishu)+loadj/2-1 ];  % 只去掉一个中间的球
       else
           loadi=[  (ceil(n*xishu)-floor(loadj/2)+1):ceil(n*xishu)+floor(loadj/2)   ];    % 去掉1个球
       end 
     end    
  end  
else
    mark1=2;              %说明收敛且没有出错！！！！
    break
end

loadj=length(loadi);   nnnum=1;
% 给定变量初值！！
x=[];r=[];
for i=1:loadj
    x=[x Fxx/1e6];   r=[r 1/1e30];
end
X2=Fxx/1e30;  Y2=(Fyy+1)/1e30;Z2=Fzz/1e30;sitay=Myy/1e8;sitaz=Mzz/1e8;       % Z2的初值随着Fxx的变化而变化，且分母也要变化，当Fxx较大时如果分母>1e8,则计算错误！！  
pppp1=[x r X2 Y2  Z2 sitay sitaz]';   %此处x r为相对球坐标系的位移！！！！！！！！！！
www3=pppp1; 
ffLOAD(datafromvb,loadi,www3);load result444;  ddd1=1e-5;    mark1=0;    % mark1要置0
end

end





% 当只有1个ball承载时，由于sita=0，sin（sita）=0，计算Y2时会出错， 所以不能考虑y2
ffLOAD(datafromvb,loadi,www3);      load result444;
if loadj==1
while( result444>0.1 )        %单次计算 循环
    ffLOAD(datafromvb,loadi,www3);      load result444;result444
    load uu; www3=uu;
    JffLOAD(datafromvb,loadi,www3);    load Jzz;
    wwwoneball=[( [www3(1:2*loadj+1);www3(2*loadj+3)]  -Jzz([1:2*loadj+1 2*loadj+3],[1:2*loadj+1 2*loadj+3])\[zzzz3(1:2*loadj+1);zzzz3(2*loadj+3)] ); ];
    www3=[wwwoneball(1:2*loadj+1);Y2;wwwoneball(2*loadj+2);sitay;sitaz];
    for i=1:loadj
        if imag(www3(i))~=0
            mark1=1;break
        end
    end
    if mark1==1
        break
    end
    if result444>1e60
            mark1=1
            break
    end
    if mark1==1
        break
    end
    nnnum=nnnum+1
    if nnnum>60
        mark1=1;    %说明该载荷下 数为 loadj 无法收敛！！！  出错!
    break
    end
end

load q1q2a1a2;
for i=1:loadj
    if   Q1(i)>(sqrt(Fxx^2+Fyy^2+Fzz^2)+Fz(i)) |  Q1(i)<1 | Q2(i)<1 | a1(i)<0 | a2(i)<0       %%判断是否出错！！（即滚动体径向位移是否为负 或者接触角是否<0 ）
        mark1=1;
        break
    end
end
if mark1==1      % 表示计算出错,此时要去掉 承载  滚动体个数!!!!
loadj=0; loadi=[];
else
    mark1=2;              %说明收敛且没有出错！！！！
end
loadj=length(loadi); 
end

save www3; save loadi;





