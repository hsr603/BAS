clc
clear all
close all

% dbstop if error


if 1  %参数设置
    global L p F_limit U_limit Area NVAR  my_rule_opst_num NIND
    D=Data10;
    Area=D.Area; p=D.p;U_limit=D.U_limit;F_limit=D.F_limit;
    T_rank=D.T; %避免跟ST中的T重复
    MAXGEN=99;
    NIND=1;NVAR=size(T_rank,2);PRECI=1;
    BaseV=ones(1,NVAR)*max(size(Area));
    rule=ones(1,NVAR);
    [A_min,A_max]=size(Area);
    my_rule_opst_num=rule;
    for i=1:size(D.Con,2) %计算杆件长度  设置成全局变量 尺寸优化 拓扑属性不变 算一次就可以
        H=D.Con(:,i);%得到第i根杆件的编号
        C=D.Coord(:,H(2))-D.Coord(:,H(1)); %杆件坐标差
        L(i)=norm(C);%norm(C)是在求这根杆件的长度
    end
    Tpr_a=1-D.Re;f=find(Tpr_a); %这一步是为了将位移转换而用到的东西
    clear C H Tpr_a
end
tic

if 1 %初始化b
    Initial_Chrom=crtbp(NIND,BaseV);
    Initial_Chrom=Initial_Chrom+1;    %+1是因为上面那个函数本身的缺陷
    Initial_Chrom; %分组后的A
    Initial_A=my10t25(Initial_Chrom,T_rank);
    x(1,:)=Initial_Chrom;
    A(1,:)=Initial_A;
end

if 1 %初始化 a)计算目标函数 
    for i=1:size(Initial_A,1)  %计算目标函数
        D.A=Initial_A(i,:)';
        [~,~,ObjV(i,:)]=my_objv(D);
    end
    clear Initial_Chrom Initial_FitnV Initial_A
end
if 1 % 一个测试模块
    x_best=x(1,:);
    y_best=ObjV(1,:);
end
d=1;
for gen=1:MAXGEN
    dir=rand(NVAR,1);  dir=dir/norm(dir);
    d0=1; %离散变量问题d0最小取值
    step=1;
    %****生成xl和xr
    d=d0+d^(gen-1)*0.95;
    xl=x(gen,:)+d*dir';xl=round(xl);
    xr=x(gen,:)-d*dir';xr=round(xr);
    %****修正
    if min(xl)<A_min
     xl(xl<A_min)=A_min;
    end
    if min(xr)<A_min
     xr(xr<A_min)=A_min;
    end
    %***下一步方向
xl(xl>A_max)=A_max;xr(xr>A_max)=A_max;
Al=my10t25(xl,T_rank);Ar=my10t25(xr,T_rank);
D.A=Al';[~,~,y_left]=my_objv(D);
D.A=Ar';[~,~,y_right]=my_objv(D);  
x(gen+1,:)=x(gen,:)+step*sign(y_left(6)-y_right(6))*dir';
 %****修正x
 x(gen+1,:)=round(x(gen+1,:));
    if min(x(gen+1,:))<A_min
     x(gen+1,x(gen+1,:)<A_min)=A_min;
    end
    if min(x)<A_min
     x(gen+1,x<A_min)=A_min;
    end
    x(gen+1,x(gen+1,:)>A_max)=A_max;
    A=my10t25(x(gen+1,:),T_rank);
    D.A=A';
    [~,~,ObjV(gen+1,:)]=my_objv(D);
    step=step*0.95;
    
    if 1 %新模块
        if ObjV(gen+1,6)<y_best(1,6)
            x_best=x(gen+1,:);
            y_best=ObjV(gen+1,:);
        end
    end
end  %开始循环优化
toc
disp('my_op.m end')




