module SBModel
    use global
    implicit none

    contains
    !陕北产流模型
    subroutine SBModel_RP(x,y,P1,E1,Output)
        implicit none
        real::x(:),y(:),Output(:,:),P1(:),E1(:)
        real::mu,KC,f0,fc,B,FB,CS,KK,xx,k,W0,W,EP
        real::PE,dW1,tt,f,fmm,R1,R2,R,DC,W1,DQT,n
        integer::i,t,nf,l
        real,allocatable::RS(:),RG(:)
        real::AA

        KC=x(1)
        fc=x(2)
        B=x(3)
        FB=0.01
        kk=x(4)

        allocate(RS(size(Q)),RG(size(Q)))
        RS=0.0;RG=0.0
        f0=y(1)+fc;W0=y(2);W=W0

        do t=2,size(Q)

            EP=E1(t)*KC;PE=P1(t)-EP
            if (PE<0.0)PE=0.0

            dW1=1.0;tt=W/f0                      
            do while(dW1>0.05)
              W1=fc*tt+1/kk*(1-exp(-kk*tt))*(f0-fc)
              dW1=abs(W1-W)
              f=f0-kk*(W1-fc*tt)
              tt=tt+abs(W1-W)/f
              if(tt>W/fc)then
                tt=W/fc
                exit
              end if
            end do
            !f=B**2*(1-(1+A*W/B**2)**0.5)/W+A

            f=f0-kk*(W-fc*tt)
            !计算地表径流
            if (abs(PE)>1E-5)then
              if(PE>f)then
                R2=PE-f
              else
                R2=0.0
              end if
            else
              R2=0.0
            end if
            W=W+P1(t)-R2-EP

            if(W<0.0)W=0.0
            if(PE>0.0)then
              R1=PE*FB
            else
              R1=0.0
            end if           
             RS(t)=R1+R2*(1.0-FB)          

      end do
      Output(:,1)=RS

    end subroutine

    !陕北汇流模型
    subroutine SBModelRouting(RS,RI,RG,x,x1,Q1,QC)
        implicit none
        real::Q1,x1
        real::RS(:),RI(:),RG(:),QC(:),x(:)
        integer::i,j,t,NN
        real::CS,n,k,xx
        real,allocatable::QS(:),R(:)

        CS=x(1)
        n=x1
        k=x(3)
        xx=x(4)
        NN=size(QC)
        allocate(QS(NN),R(NN))
        R=RS+RG+RI;QS(1)=Q1

        do t=1,NN-1
            QS(t+1)=CS*QS(t)+(1-CS)*R(t)*AA/3.6/periods
        end do

        call Muskingum(k,xx,n,QS,QC)
    end subroutine

    !分段马斯京根汇流
   subroutine Muskingum(k,xx,n1,QS,QC)
      implicit none
      real::k,xx,xxl,C0,C1,C2,n1,kl,dt
      integer::m,n,i,j,t,TT
      real::QS(:),QC(:)
      real,allocatable::CT(:),C(:)
      real,allocatable::QCP(:),QC1(:,:),QC2(:),QC3(:)

      n=floor(n1)
      kl=k
      dt=1.0
      xxl=xx
      C0=(0.5*dt-kl*xxl)/(0.5*dt+kl-kl*xxl)
      C1=(0.5*dt+kl*xxl)/(0.5*dt+kl-kl*xxl)
      C2=(-0.5*dt+kl-kl*xxl)/(0.5*dt+kl-kl*xxl) 
      TT=size(QS)
      
       !河道汇流采用马斯京根法
       allocate(QC1(TT,n+1))
       QC1=Q(1);QC1(:,1)=QS
       if(n>0)then

           allocate(QC2(TT),QC3(TT))
           QC3=0.0;QC2=QS
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
              
    end subroutine Muskingum
end module