function qiujieSPEED1(datafromvb )                  %最速下降法求解非线性方程组

load loadi;
nn34=1; ij=1; speedchange=1e-6; result1134=1e4;  min3=1e50;
for i=1:3*loadj+2
    dd3d(i)=0;
end
load q1q2a1a2;
for i=1:loadj
    Q1(i)=q1q2a1a2(i);    Q2(i)=q1q2a1a2(loadj+i);  a1(i)=q1q2a1a2(2*loadj+i);  a2(i)=q1q2a1a2(3*loadj+i);
end
for i=1:loadj
    if Q1(i)==max(Q1)
       markQ1=i;          % 找到受载最大的球
    end
end

for i=1:loadj
beitajiao(i)=atan(sin(a1(i))/(cos(a1(i))+Dw/Dm));
Wo(i)=W2*(1-Dw/Dm*cos(a2(i)))*(cos(a1(i))+tan(beitajiao(i))*sin(a1(i)))/(( 1-Dw/Dm*cos(a2(i)))*(cos(a1(i))+tan(beitajiao(i))*sin(a1(i)))+ ( 1+Dw/Dm*cos(a1(i)))*(cos(a2(i))+tan(beitajiao(i))*sin(a2(i)))    );
Wx(i)=-W2*(1-Dw/Dm*cos(a2(i)))*(1+Dw/Dm*cos(a1(i)))/(( 1-Dw/Dm*cos(a2(i)))*(cos(a1(i))+tan(beitajiao(i))*sin(a1(i)))+ ( 1+Dw/Dm*cos(a1(i)))*(cos(a2(i))+tan(beitajiao(i))*sin(a2(i)))    )/(Dw/Dm);
Wz(i)=-Wx(i)*tan(beitajiao(i));
end

if loadj==datafromvb(1)
    pppp34=[ Wo Wx Wz ]';%未知变量
else
Wononload=W2*(1-Dw/Dm)/((1-Dw/Dm)+(1+Dw/Dm));
Wxnonload=-W2*(1-Dw/Dm)*(1+Dw/Dm)/(( 1-Dw/Dm)+ ( 1+Dw/Dm) )/(Dw/Dm);
Wznonload=0;
pppp34=[ Wo Wx Wz Wononload Wxnonload]';%未知变量
end

nn34=1;  www34=pppp34;
while( result1134>0.1 )
    ffSPEED1(datafromvb,www34);    load result1134;aaa34=result1134;
    if result1134<min3
    wwmin34=www34; save wwmin34;min3=result1134
    end
    q3=www34;    
    for ij=1:loadj
    q3(ij)=www34(ij)*(1+speedchange);  ffSPEED1(datafromvb,q3);load result1134;bbb=result1134;   dd3d(ij)=(bbb-aaa34)/(speedchange*www34(ij)); q3=www34;
    q3(loadj+ij)=www34(loadj+ij)*(1+speedchange);   ffSPEED1(datafromvb,q3);load result1134;bbb=result1134;    dd3d(loadj+ij)=(bbb-aaa34)/(speedchange*www34(loadj+ij));  q3=www34;
    q3(2*loadj+ij)=www34(2*loadj+ij)*(1+speedchange);   ffSPEED1(datafromvb,q3);load result1134;bbb=result1134;   dd3d(2*loadj+ij)=(bbb-aaa34)/(speedchange*www34(2*loadj+ij));  q3=www34;
    end
    if loadj==datafromvb(1)
    else
       q3(3*loadj+1)=www34(3*loadj+1)*(1+speedchange);   ffSPEED1(datafromvb,q3);load result1134;bbb=result1134;   dd3d(3*loadj+1)=(bbb-aaa34)/(speedchange*www34(3*loadj+1));  q3=www34;
       q3(3*loadj+2)=www34(3*loadj+2)*(1+speedchange);   ffSPEED1(datafromvb,q3);load result1134;bbb=result1134;   dd3d(3*loadj+2)=(bbb-aaa34)/(speedchange*www34(3*loadj+2));  q3=www34;
    end
    
    dd=0;
    for j=1:(length(q3))
        dd=dd+dd3d(j)^2;
    end
    suanzi3=aaa34/dd;
    qqq=www34;
    for k=1:length(q3)
    qqq(k)=www34(k)-(suanzi3*dd3d(k));  
    end
    www34=qqq'
    save www34;      nn34=nn34+1
    if nn34>1              %计算10次！！
        break
    end
end
save wwmin34;







