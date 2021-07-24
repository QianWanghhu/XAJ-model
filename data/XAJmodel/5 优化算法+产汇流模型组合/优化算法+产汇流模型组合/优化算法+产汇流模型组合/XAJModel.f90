module XAJModel
    use global
    implicit none

    contains
    !新安江产流模型
    subroutine XAJ_RP(x,y,P1,E1,Output)
        implicit none
        real::x(:),y(:),Output(:,:),P1(:),E1(:)
        real,allocatable::EP(:),WL1(:),WU1(:),WD1(:),RSS(:),RII(:),RGG(:)
        real::S0                 !时段初的自由水蓄量 
        real::FR0                !时段初的产流面积比例
        real::im=0.197          !不透水面积占比
        real::WDM,WUM,WLM,SM,b,EX,kc,C,KI,KG,DWT,SMM1,SM1
        real::WL,WU,WD,WM,WMM,MS,EU,EL,ED
        real::UE,PE,W0,SS,QD
        real::Rb,R,a,AU
        real::RS,RI,RG,KIt,KGt,FR,RS1,RG1,RI1,S
        integer::t,i,j,N,periods,NN
        real::C0,C1,C2,DQT
        allocate(EP(size(Q)))

        NN=size(Q)
        allocate(WL1(NN),WU1(NN),WD1(NN),RSS(NN),RII(NN),RGG(NN))

        WDM=x(1)                 !深层土壤水含水容量
        WUM=x(2)                 !上层土壤水含水容量
        WLM=x(3)                 !下层土壤水含水容量
        SM=x(4)                  !流域平均自由水蓄积容量
        b=x(5)                   !反映流域包气带蓄水容量分布的不均匀性
        EX=x(6)                  !自由水蓄量分布曲线指数
        kc=x(7)                  !蒸发折算系数
        C=x(8)                   !蒸发扩散系数
        KI=x(9)/(x(9)+x(10))*0.8     !壤中流出流系数
        KG=x(10)/(x(9)+x(10))*0.8    !地下径流出流系数

        WM=WLM+WUM+WDM           !流域平均蓄水容量
        WMM=WM*(b+1.0)           !流域单点最大蓄水量
        MS=SM*(EX+1.0)            !流域最大自由水蓄量
        EP=kc*E1                  !计算蒸散发能力

        KIt=(1.0-(1.0-KI-KG)**(1.0/(24.0*1.0)))/(1.0+KG/KI)      !划分小段后的参数值             
        KGt=KG*KIt/KI;KI=KIt;KG=KGt                   
        WU=y(1)              !下层土壤含水量
        WL=y(2)              !深层土壤含水量
        WD=y(3)              !上层土壤含水量
        S0=y(4)
        FR0=y(5)
        if(WU>WUM)then
            WL=WL+WU-WUM;WU=WUM
        end if
        if(WL>WLM)then
            WD=WD+WL-WLM;WL=WLM
        end if          
        if(WD>WDM)WD=WDM
        if(S0>SM)then
            S0=SM
        end if

        WL1(1)=WL;WD1(1)=WD;WU1(1)=WU
        RSS(1)=0.0;RII(1)=0.0;RGG(1)=0.0

        do t=2,NN

            W0=WL+WU+WD          !t时段的初始土壤含水量
            if (W0>WM)W0=WM
            PE=P1(t)-EP(t)

            !流域蒸散发，三层蒸发模型
            if (P1(t)-EP(t)<0.0) then
                if(WU+P1(t)>EP(t))then
                    EU=EP(t);EL=0.0;ED=0.0
                    WU=WU+P1(t)-EP(t)
                else 
                    EU=WU+P1(t);WU=0.0              
                    if (WL>C*WLM) then
                        EL=(EP(t)-EU)*WL/WLM;ED=0.0;WL=WL-EL
                    else
                        if (WL>C*(EP(t)-EU)) then
                            EL=C*(EP(t)-EU);ED=0.0;WL=WL-EL
                        else
                            EL=WL;WL=0.0;ED=C*(EP(t)-EU)-WL;WD=WD-ED
                          if (WD<=0) then
                            ED=ED+WD;WD=0.0
                          end if
                        end if
                    end if
                end if
            end if
            UE=EU+ED+EL
            PE=P1(t)-UE           !扣除蒸发的净雨
            if (PE<0) PE=0.0

            !流域产流
            if (PE>0) then
                Rb=im*PE                                         !不透水面积产流
                a=WMM*(1.0-(1.0-W0/WM)**(1.0/(1.0+b)))           !W0对应的纵坐标
                if ((a+PE)<WMM) then 
                    R=PE-WM+W0+WM*((1-(PE+a)/WMM)**(b+1.0))
                else
                    R=PE-(WM-W0)
                end if
            else
                R=0.0;Rb=0.0
            end if

            if (abs(R)<1E-10) R=0.0

            if(PE>0)then
              if(WU+PE-R>WUM)then
                if(WU+WL+PE-R-WUM>WLM)then
                  WU=WUM;WL=WLM
                  WD=W0+PE-R-WU-WL
                  if(WD>WDM)WD=WDM
                else
                  WL=WU+WL+PE-R-WUM;WU=WUM                 
                end if
              else
                WU=WU+PE-R
              end if
            end if
            
            RS=0.0;RI=0.0;RG=0.0
            !分水源计算
            if (R>0.0) then
                FR=R/PE
                S=S0*FR0/FR
                QD=R/FR
                N=floor(QD/5)+1                        !将每个计算时段的入流R，按5mm分成N段
                KIt=(1.0-(1.0-KI-KG)**(1.0/(N*1.0)))/(1.0+KG/KI)      !划分小段后的参数值               
                KGt=KG*KIt/KI
                QD=QD/N

                SMM1=MS
                SM1=SMM1/(1.0+EX)
               
                SS=S;RS=0.0;RI=0.0;RG=0.0;j=0
                do i=1,N                   
                    j=j+1
                    if(S>SM1)then
                      AU=SMM1
                    else
                      AU=SMM1*(1.0-(1.0-S/SM1)**(1.0/(1.0+EX)))  !S0对应的纵坐标
                    end if

                    if (QD+AU<SMM1) then 
                        RS=FR*(QD+S-SM1+SM1*((1.0-(QD+AU)/SMM1)**(EX+1.0)))+RS
                    else
                        RS=FR*(QD+S-SM1)+RS
                    end if

                    S=i*QD-RS/FR+S
                    RI=KIt*FR*S+RI
                    RG=KGt*FR*S+RG
                    S=SS+i*QD-(RS+RI+RG)/FR
                    
                    if (S>SM1) S=SM1
                    if(S<0)S=0.0
                end do
                S0=S;FR0=FR
            else
                RS=0.0;RG=0.0;RI=0.0
                RG=S0*KG*FR0;RI=S0*KI*FR0
                S0=S0-(RG+RI)/FR0
            end if
            RS=Rb+RS*(1-im);RG=RG*(1-im);RI=RI*(1-im)
            WU1(t)=WU;WL1(t)=WL;WD1(t)=WD
            RSS(t)=RS;RII(t)=RI;RGG(t)=RG

        end do
        Output(:,1)=RSS;Output(:,2)=RII;Output(:,3)=RGG
    end subroutine

    !新安江汇流模型
    subroutine XAJRouting(RS,RI,RG,x,x1,Q1,QC)
       implicit none
       real::Q1,x1
       real::RS(:),RI(:),RG(:),QC(:),x(:)
       integer::i,j,t,NN,n,L
       real::CS,CI,CG,K,xx,dt,C0,C1,C2
       real,allocatable::QS(:),QI(:),QG(:),QT(:),QQ(:),QC1(:,:),QC2(:),QC3(:)

       CS=x(1)                 !地面径流消退系数
       CI=x(2)                 !壤中流消退系数
       CG=x(3)                 !地下径流消退系数
       L=floor(x(4))
       n=floor(x1)
       K=x(6)                  !槽蓄系数
       xx=x(7)                 !流量比重系数
       dt=K                    !演算时段

       C0=(0.5*dt-K*xx)/(0.5*dt+K-K*xx)
       C1=(0.5*dt+K*xx)/(0.5*dt+K-K*xx)
       C2=(-0.5*dt+K-K*xx)/(0.5*dt+K-K*xx)
       NN=size(QC)
       allocate(QS(NN),QI(NN),QG(NN),QT(NN),QQ(NN))
       CI=CI**(1.0/24.0)
       CG=CG**(1.0/24.0)
       QS=0.0;QI=0.0;QG=0.0;QG(1)=Q(1);QT=0.0
       QQ=Q(1)
       QQ=0.0
       QG=0.0

       !地表、壤中和地下径流以线性水库汇流
       do t=2,NN
           QS(t)=RS(t)*AA/3.6/periods
           QI(t)=CI*QI(t-1)+(1.0-CI)*RI(t)*AA/3.6/periods
           QG(t)=CG*QG(t-1)+(1.0-CG)*RG(t)*AA/3.6/periods
       end do
       QT=QS+QI+QG

       !河网汇流采用线性水库
       do t=1,NN-1
          if(t-L>1)then
            QQ(t+1)=CS*QQ(t)+(1.0-CS)*QT(t-L)  
          end if            
       end do

       !河道汇流采用马斯京根法
       allocate(QC1(TT,n+1))
       QC1=Q(1);QC1(:,1)=QQ
       if(n>0)then
           allocate(QC2(TT),QC3(TT))
           QC3=Q(1);QC2=QQ
           do i=1,n
                do t=1,TT-1
                    QC3(t+1)=C0*QC2(t+1)+C1*QC2(t)+C2*QC3(t)
                end do
                QC2=QC3
           end do
           QC=QC2
       else
           QC=QS
       end if    

    end subroutine

end module