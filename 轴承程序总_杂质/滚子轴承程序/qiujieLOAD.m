function qiujieLOAD(datafromvb)                %牛顿 法求解非线性方程组 得到承载区载荷分布 非歪斜载荷求解

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
 %================== 是否允许自动删除承载滚子 ==================
% 0：调试模式，不自动删滚子，用来观察游隙变化后候选承载区是否变化
% 1：正式模式，允许删除弱承载滚子，但改为删除 Q2 最小的滚子，而不是机械删除中间滚子
allowDeleteRoller = 1;

% 内圈接触力小于这个值时，认为该滚子不稳定承载
Q2_min_threshold = 1;     % N，原程序等效用的是 1 N
residual_accept_threshold = 1e5;  % 载荷方程残差低于该值且载荷为正时，接受当前承载区
% =============================================================
% 分配承载滚子的分布
ff1(datafromvb);% 初步求解载荷分布！！
load loadn
loadj=2*loadn+3; % 为了计算正确，将静力学计算的承载滚子数增大2个！！！！
mark1=0; % mark1为计算标记，=1代表出错，=2收敛，=3退出循环1   % Attention!!! Myy不可以为零！！否则 雅克比出错！！！
if Fyy==0    
    if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[1:loadj/2-1  n-loadj/2:n];           % 只去掉一个中间的球
       else
           loadi=[1:floor(loadj/2)  n-floor(loadj/2):n];  % 去掉1个球
    end
else              
     lipianjiao=atan(Fyy/Fzz);   % Fyy 与 Fzz 的合力与Z 轴的夹角！！！ 
     xishu=lipianjiao/(2*pi);    % 当存在Fyy 时，承载滚动体个数 就不能从最下端开始了，而应该沿 合力的反方向！！！
   if (ceil(n*xishu)-(loadj/2-1))<1        % 当可以从1号球开始时
      if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[1:ceil(n*xishu)+loadj/2-1  ceil(n*xishu)+loadj/2:loadj];   % 只去掉一个中间的球
       else
           loadi=[1:ceil(n*xishu)+floor(loadj/2) ceil(n*xishu)+floor(loadj/2)+1:loadj];    % 去掉1个球
       end
     else         % 当 不 可以从1号球开始时
        if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
           loadi=[ (ceil(n*xishu)-(loadj/2-1)):ceil(n*xishu)+loadj/2 ];  % 
       else
           loadi=[  (ceil(n*xishu)-floor(loadj/2)+1):ceil(n*xishu)+floor(loadj/2)+1   ];    % 去掉1个球
       end 
     end  
end

%初值给定
pp=[];pp(2*loadj+1 )=Fzz/1e6; pp(2*loadj+2)=Fyy/1e6;pp(2*loadj+3)=Myy/1e6;pp(2*loadj+4)=0;
for i=1:loadj
    sita(i)=2*pi*(loadi(i))/n; 
   if (pp(2*loadj+1 )*cos(sita(i))+pp(2*loadj+2)*sin(sita(i))-deltar0/2)/2>0
       pp(i)=(pp(2*loadj+1 )*cos(sita(i))+pp(2*loadj+2)*sin(sita(i))-deltar0/2)/2;   %保证初值给定在正确的范围！！保证delta1 delta2都为正值！！
    else
        pp(i)=10e-6;
    end
    pp(loadj+i )=0;
end

% 求解
www1=pp';result111=1e4; minresult111=1e1000;  mm=1; temp11=[];
while(loadj>1 & mark1~=2 )           % 循环1   loadj 承载滚子个数循环计算
    
while( result111>Myy^2)                % 循环2   不包括轴承倾覆力矩的计算！！！ 作为循环3的初始解！！
    ff2(datafromvb,loadi,www1);    load zz1;  load result111;result111
    Jff2(datafromvb,loadi,www1);
    load Jzz;
    Jzz1=Jzz([1: 2*loadj+2],[1: 2*loadj+2]);       % 分开求解，！！！！
    zz11=zz1([1: 2*loadj+2]); zz11=zz11(:);
    temp11=Jzz1\zz11;
    if any(~isfinite(temp11))
        mark1=1;
        break
    end
    oldwww1=www1(:); oldresult111=result111; acceptedStep=0; alpha=1;
    while alpha>=1/32
        wwww=oldwww1(1:(2*loadj+2)); wwww=wwww(:)-alpha*temp11(:);
        www1=[wwww; oldwww1((2*loadj+3):(2*loadj+4)); ];
        ff2(datafromvb,loadi,www1); load result111;
        if isfinite(result111) && result111<oldresult111*1.2 && all(isfinite(www1))
            acceptedStep=1;
            break
        end
        alpha=alpha/2;
    end
    if acceptedStep==0
        www1=oldwww1; result111=oldresult111; mark1=1;
        break
    end
    save www1
    mm=mm+1
    if mm>50
    break
    end
end


nn=1;
while( result111>0.1 )         % 循环3     包括轴承倾覆力矩的计算！！！  
    ff2(datafromvb,loadi,www1);   
     %load Q2;
    % for i=1:loadj                %  当循环3收敛完成时，还要判断是否出错（与内圈接触力矩为0或负值的话出错！！）！！！
        %if Q2(i)<=1
            %mark1=1; break
        %end
     %end
    % if mark1==1
       %  break;
    % end
    load Q2;
weakIndex = find(Q2 <= Q2_min_threshold);

if ~isempty(weakIndex)
    fprintf('\n[警告] 当前候选承载区存在弱承载滚子：\n');
    fprintf('loadj = %d\n', loadj);
    fprintf('loadi = ');
    disp(loadi);
    fprintf('Q2 = ');
    disp(Q2);
    if allowDeleteRoller == 1
        mark1 = 1;
        break;
    else
        fprintf('[调试模式] 检测到弱承载滚子，但暂不删除，用于观察候选承载区。\n');
    end
end
    load zz1;    yyy=zz1; load result111;result111
     if nn>50                            % 如果计算了很久还没有 收敛的迹象 ，说明成载滚动体数目不正确，退出！！！  
        if result111<minresult111
           minresult111=result111
        end
        if result111/minresult111>1e3
                 mark1=1
                 break    
        end
        if mark1==1
            break;
        end
     end
     if nn>100                            % 如果计算了很久还没有 收敛的迹象 ，说明成载滚动体数目不正确，退出！！！  
        if result111>1e5
                 mark1=1
                 break    
        end
        if mark1==1
            break;
        end
     end
    Jff2(datafromvb,loadi,www1);
    load Jzz;
    
    Jzz1=Jzz;       % 分开求解，若整体致使数量级发生错误， 无法收敛！！！！
    zz11=zz1; zz11=zz11(:);
    temp12=Jzz1\zz11;
    if any(~isfinite(temp12))
        mark1=1;
        break
    end
    oldwww1=www1(:); oldresult111=result111; acceptedStep=0; alpha=1;
    while alpha>=1/32
        wwww=oldwww1(1:(2*loadj+3)); wwww=wwww(:)-alpha*temp12(:);
        www1=[wwww;oldwww1(2*loadj+4)];
        ff2(datafromvb,loadi,www1); load result111;
        if isfinite(result111) && result111<oldresult111*1.2 && all(isfinite(www1))
            acceptedStep=1;
            break
        end
        alpha=alpha/2;
    end
    if acceptedStep==0
        www1=oldwww1; result111=oldresult111; mark1=1;
        break
    end
    save www1
    nn=nn+1
    if nn>150
        ff2(datafromvb,loadi,www1)
    break
    end
end

%load Q1;load Q2;
   %for i=1:loadj                %  当循环3收敛完成时，还要判断是否出错（可能 接触力 为负值或倾角与力矩反号！！）！！！
        %if Q2(i)<=1
            %mark1=1; break
        %if Myy>0.1
           % if www1(loadj+i)<0
             %   mark1=1; break
         %   end
          %  if www1(2*loadj+3)<0
            %    mark1=1; break
        %    end
     %   else
     %   end             
 %  end

load Q1; load Q2;

weakIndex = find(Q2 <= Q2_min_threshold);

if ~isempty(weakIndex)
    fprintf('\n[最终检查] 以下滚子内圈接触力过小：\n');
    fprintf('滚子编号 loadi = ');
    disp(loadi(weakIndex));
    fprintf('对应 Q2 = ');
    disp(Q2(weakIndex));

    if allowDeleteRoller == 1
        mark1 = 1;
    else
        fprintf('[调试模式] 最终检查发现弱承载滚子，但不删除。\n');
    end
end

if Myy > 0.1
    for i = 1:loadj
        if www1(loadj+i) < 0
            mark1 = 1; 
            break
        end
        if www1(2*loadj+3) < 0
            mark1 = 1; 
            break
        end
    end
end
if result111<=1e-1 & mark1~=1
   mark1=2;       % 代表没有出错且已经收敛！！  
elseif mark1~=1 && isfinite(result111) && result111<=residual_accept_threshold && all(isfinite(Q2)) && min(Q2)>Q2_min_threshold
   fprintf('\n[收敛提示] 残差 %.6g 高于严格阈值但载荷为正，接受当前承载区。\n', result111);
   mark1=2;
else


%if Fyy==0    
   % if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
          % loadi=[loadi(1:loadj/2-1) loadi(loadj/2+1:loadj)];           % 只去掉一个中间的球
    %   else
          % loadi=[loadi(1:floor(loadj/2))  loadi(floor(loadj/2)+2:loadj)];  % 去掉1个球
   % end
%else              
    % lipianjiao=atan(Fyy/Fzz);   % Fyy 与 Fzz 的合力与Z 轴的夹角！！！ 
    % xishu=lipianjiao/(2*pi);    % 当存在Fyy 时， 去掉承载滚动体个数时 就不能从最下端开始了，而应该沿 合力的反方向！！！
     %if (ceil(n*xishu)-(loadj/2-1))<1        % 当可以从1号球开始时
    %  if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
        %   loadi=[loadi(1:ceil(n*xishu))  loadi(ceil(n*xishu)+1:ceil(n*xishu)+loadj/2-1)  loadi(ceil(n*xishu)+loadj/2+1:loadj)];   % 只去掉一个中间的球
      % else
          % loadi=[loadi(1:ceil(n*xishu))  loadi(ceil(n*xishu)+1:ceil(n*xishu)+floor(loadj/2)) loadi(ceil(n*xishu)+floor(loadj/2)+2:loadj)];    % 去掉1个球
       %end
   %  else         % 当 不 可以从1号球开始时
     %   if mod(loadj/2,1)==0        %说明滚动体个数为偶数！！
     %      loadi=[ (ceil(n*xishu)-(loadj/2-1)):ceil(n*xishu)+loadj/2-1 ];  % 只去掉一个中间的球
     %  else
     %      loadi=[  (ceil(n*xishu)-floor(loadj/2)+1):ceil(n*xishu)+floor(loadj/2)   ];    % 去掉1个球
      % end 
   %  end   
%end 
%loadi
%loadj=length(loadi);  

% ============= 自动删滚子逻辑：改为删除 Q2 最小的那个滚子 =============
if allowDeleteRoller == 1

    if exist('Q2','var') == 0
        load Q2;
    end

    if length(Q2) == length(loadi)
        [q2min, idxDelete] = min(Q2);
    else
        % 如果 Q2 没有正常更新，则默认删除候选区边缘滚子
        idxDelete = 1;
        q2min = NaN;
    end

    fprintf('\n[自动调整承载滚子数]\n');
    fprintf('删除前 loadj = %d\n', loadj);
    fprintf('删除前 loadi = ');
    disp(loadi);
    fprintf('删除滚子编号 = %d, Q2 = %.6f N, result111 = %.6g\n', loadi(idxDelete), q2min, result111);

    loadi(idxDelete) = [];
    loadj = length(loadi);

    fprintf('删除后 loadj = %d\n', loadj);
    fprintf('删除后 loadi = ');
    disp(loadi);

else
    fprintf('\n[调试模式] 当前不允许自动删滚子，保留候选承载区。\n');
    fprintf('当前 loadj = %d\n', loadj);
    fprintf('当前 loadi = ');
    disp(loadi);
    mark1 = 2;
    break;
end
% =====================================================================

  pp=[];
  pp(2*loadj+1 )=Fzz/1e6; pp(2*loadj+2)=Fyy/1e6;pp(2*loadj+3)=Myy/1e6;pp(2*loadj+4)=0;
  for i=1:loadj
    sita(i)=2*pi*(loadi(i))/n; 
    if (pp(2*loadj+1 )*cos(sita(i))+pp(2*loadj+2)*sin(sita(i))-deltar0/2)/2>0
       pp(i)=(pp(2*loadj+1 )*cos(sita(i))+pp(2*loadj+2)*sin(sita(i))-deltar0/2)/2;   %保证初值给定在正确的范围！！保证delta1 delta2都为正值！！
    else
        pp(i)=10e-6;
    end
    pp(loadj+i )=0;
  end
  www1=pp';  result111=1e4;   mm=1; mark1=0;  
end
end

% 当只有1个滚子承载时，由于sita=0，sin（sita）=0，计算Y2时会出错， 所以不能考虑y2，程序要重新改变

if loadj==1
     pp=[];
  pp(2*loadj+1 )=Fzz/1e4; pp(2*loadj+2)=Fyy/1e6;pp(2*loadj+3)=Myy/1e6;pp(2*loadj+4)=0;
  for i=1:loadj
    sita(i)=2*pi*(loadi(i))/n; 
    pp(i)=(pp(2*loadj+1 )*cos(sita(i))+pp(2*loadj+2)*sin(sita(i))-deltar0/2)/2;
    pp(loadj+i )=0;
  end
  www1=pp';  result111=1e4;   mm=1; mark1=0;  
    
while( result111>Myy^2)                % 循环2   不包括轴承倾覆力矩的计算！！！ 作为循环3的初始解！！
    
    ff2(datafromvb,loadi,www1);    load zz1;  load result111;result111

    Jff2(datafromvb,loadi,www1);
    load Jzz;
    Jzz1=Jzz([1: 2*loadj+1],[1: 2*loadj+1]);       % 分开求解，！！！！
    zz11=zz1([1: 2*loadj+1]); zz11=zz11(:);
    temp11=Jzz1\zz11;
    for i=1:loadj+1
        if temp11(i)~=0
           mark1=3; break    
        end
    end
    if mark1~=3        % 即如果 temp11 全为0(此时收敛不起作用了) ，则退出循环2
        break
    end
    
    wwww=www1([1:2*loadj+1 ])-temp11;
    www1=[wwww; www1(2*loadj+2:2*loadj+4); ]
    save www1
    mm=mm+1
    if mm>50
    break
    end
end

end

save www1; save loadi;

% 非承载区载荷计算：认为非承载区受力完全相同！ Q2=0;

