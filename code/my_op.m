clc
clear all
close all

% dbstop if error


if 1  %��������
    global L p F_limit U_limit Area NVAR  my_rule_opst_num NIND
    D=Data10;
    Area=D.Area; p=D.p;U_limit=D.U_limit;F_limit=D.F_limit;
    T_rank=D.T; %�����ST�е�T�ظ�
    MAXGEN=99;
    NIND=1;NVAR=size(T_rank,2);PRECI=1;
    BaseV=ones(1,NVAR)*max(size(Area));
    rule=ones(1,NVAR);
    [A_min,A_max]=size(Area);
    my_rule_opst_num=rule;
    for i=1:size(D.Con,2) %����˼�����  ���ó�ȫ�ֱ��� �ߴ��Ż� �������Բ��� ��һ�ξͿ���
        H=D.Con(:,i);%�õ���i���˼��ı��
        C=D.Coord(:,H(2))-D.Coord(:,H(1)); %�˼������
        L(i)=norm(C);%norm(C)����������˼��ĳ���
    end
    Tpr_a=1-D.Re;f=find(Tpr_a); %��һ����Ϊ�˽�λ��ת�����õ��Ķ���
    clear C H Tpr_a
end
tic

if 1 %��ʼ��b
    Initial_Chrom=crtbp(NIND,BaseV);
    Initial_Chrom=Initial_Chrom+1;    %+1����Ϊ�����Ǹ����������ȱ��
    Initial_Chrom; %������A
    Initial_A=my10t25(Initial_Chrom,T_rank);
    x(1,:)=Initial_Chrom;
    A(1,:)=Initial_A;
end

if 1 %��ʼ�� a)����Ŀ�꺯�� 
    for i=1:size(Initial_A,1)  %����Ŀ�꺯��
        D.A=Initial_A(i,:)';
        [~,~,ObjV(i,:)]=my_objv(D);
    end
    clear Initial_Chrom Initial_FitnV Initial_A
end
if 1 % һ������ģ��
    x_best=x(1,:);
    y_best=ObjV(1,:);
end
d=1;
for gen=1:MAXGEN
    dir=rand(NVAR,1);  dir=dir/norm(dir);
    d0=1; %��ɢ��������d0��Сȡֵ
    step=1;
    %****����xl��xr
    d=d0+d^(gen-1)*0.95;
    xl=x(gen,:)+d*dir';xl=round(xl);
    xr=x(gen,:)-d*dir';xr=round(xr);
    %****����
    if min(xl)<A_min
     xl(xl<A_min)=A_min;
    end
    if min(xr)<A_min
     xr(xr<A_min)=A_min;
    end
    %***��һ������
xl(xl>A_max)=A_max;xr(xr>A_max)=A_max;
Al=my10t25(xl,T_rank);Ar=my10t25(xr,T_rank);
D.A=Al';[~,~,y_left]=my_objv(D);
D.A=Ar';[~,~,y_right]=my_objv(D);  
x(gen+1,:)=x(gen,:)+step*sign(y_left(6)-y_right(6))*dir';
 %****����x
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
    
    if 1 %��ģ��
        if ObjV(gen+1,6)<y_best(1,6)
            x_best=x(gen+1,:);
            y_best=ObjV(gen+1,:);
        end
    end
end  %��ʼѭ���Ż�
toc
disp('my_op.m end')




