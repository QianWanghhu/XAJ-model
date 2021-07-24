module HBModel
    use global
    implicit none

    contains
    !河北产流模型
    subroutine HBModel_RP(x,y,P1,E1,Output)
      implicit none
      real::x(:),y(:),Output(:,:),P1(:),E1(:)
      real::FM,N,Fc,WM,b,F,U,WMM,W0,DC,A,DQT,W,KC,IM,I0,I1,X1,K,PA0,WU,WL,WUM,WLM,EU,EL,WMX
      integer::NN
      integer::t,i,TG,TS
      real,allocatable::PE(:),R(:),RG(:)
      real::s1,s2,s,QS1,QG1,B1,CC,ss,ss1,ss2,FA

      U=x(1)
      FC=x(2)
      FM=x(3)
      N=x(4)
      WM=x(5)
      WMX=x(6)
      B=x(7)
      X1=x(8)
      NN=size(Q)

      allocate(R(NN),RG(NN),PE(NN))
      
      WMM=WM*(B+1);WUM=WM*WMX;WLM=WM-WUM
      IM=X1*WM;R=0.0;RG=0.0
      I0=y(1);PA0=y(2)      
      W0=PA0;WU=PA0*WMX;WL=PA0-WU

      !计算产流量
      do t=2,NN
           
          F=(P1(t)-P1(t)**(1.0+N)/FM**N/(1.0+N))*exp(-u*I0)+FC
          FA=P1(t)-P1(t)**(1.0+N)/FM**N/(1.0+N)
          
          if(P1(t)>0.0)then
              F=FC+1.0;ss1=1.0
              do while(ss1>1E-3)
                  ss2=FA*exp(-u*I0)*(1-exp(-u*F))/(u*F)+FC
                  ss1=abs(ss2-F)
                  F=ss2
              end do
           end if
                               
          !地表产流
          R(t)=P1(t)-F    
          if(R(t)<0) R(t)=0
          if(P1(t)<F)then
              F=P1(t)
          else
              F=F
          end if

          !地下产流
          if(F>0.0)then
              if(W0<WM)then
                  A=WMM*(1.0-(1.0-W0/WM)**(1.0/(1.0+B)))
              else
                  A=WMM
              end if
              if((F+A)>WMM)then
                  RG(t)=F-(WM-W0)
              else
                  RG(t)=F-(WM-W0)+WM*(1.0-(F+A)/WMM)**(1.0+b)
              end if
          else
              RG(t)=0.0
          end if


          !计算蒸散发及土湿
          if(WU+F-E1(t)>0)then
              EU=E1(t);EL=0.0;WU=WU+F-E1(t)
          else
              EU=WU+F;WU=0.0;EL=(E1(t)-EU)*WL/WLM;WL=WL-EL
          end if

          if(F>0)then
              if(WU+F-RG(t)>WUM)then
                  WL=WL+WU+F-RG(t)-WUM;WU=WUM
                  if(WL>WLM)WL=WLM
              else
                  WU=WU+F-RG(t)
              end if
          end if

          W0=WL+WU
          !时段末状态量
          if(W0<0)W0=0
          K=1-E(t)/IM;I1=K*(I0+P1(t)-R(t))
          if(I1<0)I1=0
        end do

        Output(:,1)=R;Output(:,3)=RG
    end subroutine

    !河北汇流模型
    subroutine HBModelRouting(RS,RI,RG,x,x1,Q1,QC)
        implicit none
        real::Q1,x1
        real::RS(:),RI(:),RG(:),QC(:),x(:)
        integer::i,j,t,NN,n
        real::AS,AG,BS,BG,ss,ss1,B1,ss2
        integer::TS,TG
        real,allocatable::QS(:),QG(:)

        AS=x(1)
        AG=x(2)
        BS=x(3)
        BG=x(4)     
        TS=floor(x(5))
        TG=floor(x(6))
        NN=size(QC)
        AS=AS*AA**1.2/1.0/60.0
        AG=AG*AA**1.2/1.0/60.0
        allocate(QS(NN),QG(NN))
        QS=0.0;QG=Q1;QC=Q1
        RG=RG+RI

        !坡面、河网和河道统一采用非线性汇流模型
        do t=1,NN-1

           !地表汇流计算      
            if(t<TS+1)then
                QS(t)=0.0
            else  
                if(RS(t-TS+1)>1E-4) then
                  B1=AS/(1.0-BS)
                  ss=RS(t-TS+1)*AA/3.6/1.0-QS(t)/2.0+B1*QS(t)**(1.0-BS)
                  QS(t+1)=QS(t)+1
                  ss1=abs(QS(t+1)-QS(t))
                  do while(ss1>1E-3)
                    ss2=ss/(B1*QS(t+1)**(-BS)+0.5)
                    ss1=abs(ss2-QS(t+1))
                    QS(t+1)=ss2
                  end do                         
              else             
                  QS(t+1)=QS(t)*(1+BS/AS*QS(t)**BS)**(-1/BS)
              end if
            end if

            !地下汇流计算
            if(t<TG+1)then 
                QG(t+1)=QG(t)*(1+BG/AG*QG(t)**BG)**(-1/BG)
            else
                if(RG(t-TG+1)>1E-4)then
                  B1=AG/(1.0-BG)
                  ss=RG(t-TG+1)*AA/3.6/1.0-QG(t)/2.0+B1*QG(t)**(1.0-BG)
                  QG(t+1)=QG(t)+1
                  ss1=abs(QG(t+1)-QG(t))
                  do while(ss1>1E-3)
                    ss2=ss/(B1*QG(t+1)**(-BG)+0.5)
                    ss1=abs(ss2-QG(t+1))
                    QG(t+1)=ss2
                  end do                            
                else
                  QG(t+1)=QG(t)*(1+BG/AG*QG(t)**BG)**(-1/BG)
              end if 
            end if                  
            if(QS(t+1)<0)QS(t+1)=0.0
            if(QG(t+1)<0)QG(t+1)=0.0
        end do

        QC=QS+QG      
    end subroutine

end module