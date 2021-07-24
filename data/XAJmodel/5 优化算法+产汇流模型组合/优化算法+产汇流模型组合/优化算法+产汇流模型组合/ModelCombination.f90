module ModelCombination
    use global
    use XAJModel
    use HBModel
    use SBModel
    implicit none

    contains
    subroutine ParameterScale(n)
        implicit none
        integer::n,i,j

        select case(PRModel)
        case(1)
            NRP=10
        case(2)
            NRP=8
        case(3)
            NRP=4
        end select

        select case(RoutingModel)
        case(1)
            NRR=7
        case(2)
            NRR=6
        case(3)
            NRR=4
        end select

        allocate(UpB(NRR+NRP+3*IN+NP),DownB(NRR+NRP+3*IN+NP),Xbest(NRR+NRP+3*IN+NP))
        UpB=0.0;DownB=0.0;Xbest=0.0

        !Key2为1时需要输入参数范围，为1时直接输入参数数值
        if(Key2==1)then
            
            if(PRModel==1)then
                filename=trim(address)//'\需设置文件\产流模型参数\新安江模型范围.txt'
            else if(PRModel==2)then
                filename=trim(address)//'\需设置文件\产流模型参数\河北模型范围.txt'
            else
                filename=trim(address)//'\需设置文件\产流模型参数\陕北模型范围.txt'
            end if   
            open(10,file=filename)
            do i=1,NRP
                read(10,*)UpB(i),DownB(i)
            end do
            close(10)

            if(RoutingModel==1)then
                filename=trim(address)//'\需设置文件\汇流模型参数\新安江模型范围.txt'
            else if(RoutingModel==2)then
                filename=trim(address)//'\需设置文件\汇流模型参数\河北模型范围.txt'
            else
                filename=trim(address)//'\需设置文件\汇流模型参数\陕北模型范围.txt'
            end if

            open(10,file=filename)
            do i=1,NRR
                read(10,*)UpB(NRP+i),DownB(NRP+i)
            end do
            close(10)

            if(IN>0)then
                do i=1,IN
                    filename=trim(address)//'\需设置文件\上游控制站\马斯京根参数范围.txt'
                    open(10,file=filename)
                    do j=1,3
                        read(10,*)UpB(NRR+NRP+(i-1)*3+j),DownB(NRR+NRP+(i-1)*3+j)
                    end do
                    close(10)
                end do
            end if

            filename=trim(address)//'\需设置文件\汇流模型参数\马斯京根分段数范围.txt'
            open(10,file=filename)
            do j=1,NP
                read(10,*)UpB(NRR+NRP+3*IN+j),DownB(NRR+NRP+3*IN+j)
            end do
            close(10)

        else

            if(PRModel==1)then
                filename=trim(address)//'\需设置文件\产流模型参数\新安江模型.txt'
            else if(PRModel==2)then
                filename=trim(address)//'\需设置文件\产流模型参数\河北模型.txt'
            else
                filename=trim(address)//'\需设置文件\产流模型参数\陕北模型.txt'
            end if   
            open(10,file=filename)
            do i=1,NRP
                read(10,*)Xbest(i)
            end do
            close(10)

            if(RoutingModel==1)then
                filename=trim(address)//'\需设置文件\汇流模型参数\新安江模型.txt'
            else if(RoutingModel==2)then
                filename=trim(address)//'\需设置文件\汇流模型参数\河北模型.txt'
            else
                filename=trim(address)//'\需设置文件\汇流模型参数\陕北模型.txt'
            end if
            open(10,file=filename)
            do i=1,NRR
                read(10,*)Xbest(NRP+i)
            end do
            close(10)

            if(IN>0)then
                filename=trim(address)//'\需设置文件\上游控制站\马斯京根参数.txt'
                open(10,file=filename)
                read(10,*)
                do i=1,3*IN
                    if((MOD(i,3)==1).AND.(i/=3*IN))then
                        read(10,*)
                    end if
                    read(10,*)Xbest(NRR+NRP+i)
                end do
                close(10)
            end if      
            filename=trim(address)//'\需设置文件\汇流模型参数\子流域马斯京根分段数.txt'
            open(10,file=filename)
            do j=1,NP
                read(10,*)Xbest(NRR+NRP+3*IN+i)
            end do
            close(10)                             
        end if

    end subroutine

    subroutine StateBoundGet(Xbest,UpState,DownState)
        implicit none
        integer::i
        real,allocatable::Xbest(:),UpState(:),DownState(:)

        Select case(PRModel)
        case(1)
            allocate(UpState(5),DownState(5))
            UpState=(/Xbest(2),Xbest(3),Xbest(1),Xbest(4),1.0/)
            DownState=(/1E-5,1E-5,1E-5,1E-5,1E-5/)
        case(2)
            allocate(UpState(2),DownState(2))
            UpState=(/Xbest(5),Xbest(5)/)
            DownState=(/1E-5,1E-5/)
        case(3)
            allocate(UpState(2),DownState(2))
            UpState=(/100.0,500.0/)
            DownState=(/0.1,0.1/)
        end select


    end subroutine


    subroutine CombinedModel(x,y1,QC,OFV)
        
        !不同模型组合产汇流
        implicit none
        integer::i,j,N1
        real::x(:),y1(:),QC(:)
        real::OFV,DC,DWT,DQT
        real,allocatable::RunoffProduction(:,:),RunoffProductionAll(:,:),y(:),k(:),xx(:),nn(:),QM1(:),QM2(:),QC1(:)
        real,allocatable::RS(:),RI(:),RG(:),maxQ(:),maxQC(:),DW(:)

        allocate(RunoffProduction(TT,3),RunoffProductionAll(TT,3),maxQ(NF),maxQC(NF),DW(NF),RS(TT),RI(TT),RG(TT),k(IN),xx(IN),nn(IN),QM1(TT),QM2(TT),QC1(TT))
        RunoffProduction=0.0;RunoffProductionAll=0.0;QM1=0.0;QM2=0.0;QC1=0.0;QC=0.0
        if(PRModel==1)then
            N1=5
        else
            N1=2
        end if
        
        if(IN>0)then
            do i=1,IN
                k(i)=x(NRP+NRR+(i-1)*3+1);xx(i)=x(NRP+NRR+(i-1)*3+2);nn(i)=x(NRP+NRR+(i-1)*3+3)
            end do
        end if

        allocate(y(N1))
        do i=1,NF
            
            !产流模型         
            do j=1,NP
                y(:)=y1((i-1)*NP*N1+(j-1)*N1+1:(i-1)*NP*N1+j*N1);RunoffProduction=0.0;QC1=0.0
                select case(PRModel)
                case(1)
                    call XAJ_RP(x(1:NRP),y(:),P(FloodTime(i)+1:FloodTime(i+1),j),E(FloodTime(i)+1:FloodTime(i+1)),RunoffProduction(FloodTime(i)+1:FloodTime(i+1),:))
                case(2)
                    call HBModel_RP(x(1:NRP),y(:),P(FloodTime(i)+1:FloodTime(i+1),j),E(FloodTime(i)+1:FloodTime(i+1)),RunoffProduction(FloodTime(i)+1:FloodTime(i+1),:))
                case(3)
                    call SBModel_RP(x(1:NRP),y(:),P(FloodTime(i)+1:FloodTime(i+1),j),E(FloodTime(i)+1:FloodTime(i+1)),RunoffProduction(FloodTime(i)+1:FloodTime(i+1),:))
                end select
                RS=RunoffProduction(FloodTime(i)+1:FloodTime(i+1),1);RI=RunoffProduction(FloodTime(i)+1:FloodTime(i+1),2);RG=RunoffProduction(FloodTime(i)+1:FloodTime(i+1),3)

                !汇流模型 
                select case(RoutingModel)
                case(1)
                    call XAJRouting(RS(FloodTime(i)+1:FloodTime(i+1)),RI(FloodTime(i)+1:FloodTime(i+1)),RG(FloodTime(i)+1:FloodTime(i+1)),x(NRP+1:NRP+NRR),x(NRR+NRP+3*IN+j),Q(FloodTime(i)+1),QC1(FloodTime(i)+1:FloodTime(i+1)))
                case(2)
                    call HBModelRouting(RS(FloodTime(i)+1:FloodTime(i+1)),RI(FloodTime(i)+1:FloodTime(i+1)),RG(FloodTime(i)+1:FloodTime(i+1)),x(NRP+1:NRP+NRR),x(NRR+NRP+3*IN+j),Q(FloodTime(i)+1),QC1(FloodTime(i)+1:FloodTime(i+1)))     
                case(3)
                    call SBModelRouting(RS(FloodTime(i)+1:FloodTime(i+1)),RI(FloodTime(i)+1:FloodTime(i+1)),RG(FloodTime(i)+1:FloodTime(i+1)),x(NRP+1:NRP+NRR),x(NRR+NRP+3*IN+j),Q(FloodTime(i)+1),QC1(FloodTime(i)+1:FloodTime(i+1))) 
                end select

                QC(FloodTime(i)+1:FloodTime(i+1))=QC(FloodTime(i)+1:FloodTime(i+1))+QC1(FloodTime(i)+1:FloodTime(i+1))*weight(j)
            end do

            if(IN>0)then
                do j=1,IN
                    QM1(FloodTime(i)+1:FloodTime(i+1))=0.0
                    call Muskingum_Routing(QIN(:,j),k(j),xx(j),nn(j),1.0,QM1(FloodTime(i)+1:FloodTime(i+1)))
                    QM2(FloodTime(i)+1:FloodTime(i+1))=QM2(FloodTime(i)+1:FloodTime(i+1))+QM1(FloodTime(i)+1:FloodTime(i+1))
                end do
            end if

        end do

        QC=QC+QM2

        !模拟结果的目标函数计算
        do i=1,NF
            maxQ(i)=maxval(Q(FloodTime(i)+1:FloodTime(i+1)))
            maxQC(i)=maxval(QC(FloodTime(i)+1:FloodTime(i+1)))
            DW(i)=abs(sum(QC(FloodTime(i)+1:FloodTime(i+1)))-sum(Q(FloodTime(i)+1:FloodTime(i+1))))/sum(Q(FloodTime(i)+1:FloodTime(i+1)))
        end do
        DC=1-(sum((Q-QC)**2.0))/(sum((Q-sum(Q)/size(Q))**2.0))
        DQT=sum(abs(maxQ-maxQC)/maxQ)/NF
        DWT=sum(DW)/NF
        OFV=CC(1)*DQT+CC(2)*(1-DC)+CC(3)*DWT

    end subroutine

    subroutine Muskingum_Routing(QIN,k,xx,NN,dt,QC)
       implicit none
       real::k,xx,NN,dt,t
       real::QIN(:),QC(:)
       Integer::i,j,N1,TT
       real,allocatable::QC1(:,:)
       real::C0,C1,C2

       N1=ceiling(NN)
       C0=(0.5*dt-k*xx)/(0.5*dt+k-k*xx)
       C1=(0.5*dt+k*xx)/(0.5*dt+k-k*xx)
       C2=(-0.5*dt+k-k*xx)/(0.5*dt+k-k*xx)
       TT=size(QC)

       !河道汇流采用马斯京根法
       allocate(QC1(TT,N1+1))
       QC1=QIN(1);QC1(:,1)=QIN
       do i=1,N1
          do t=2,TT
              QC1(t,i+1)=C0*QC1(t,i)+C1*QC1(t-1,i+1)+C2*QC1(t-1,i)
          end do
       end do
       QC=QC1(:,N1+1)

    end subroutine

end module