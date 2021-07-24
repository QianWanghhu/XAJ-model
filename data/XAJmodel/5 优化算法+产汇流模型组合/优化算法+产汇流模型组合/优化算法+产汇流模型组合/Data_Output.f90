module Data_Output
    use global
    implicit none

    contains
    !结果输出
    subroutine DataOutput(QC)
        implicit none
        character(len=2000)::Number,Cline
        real::QC(:)
        integer::i,j,N1
        
        !输出最优产流参数
        filename=trim(address)//'\结果\最优产流参数.txt'
        open(10,file=filename)
        do i=1,NRP
            write(10,*)Xbest(i)
        end do
        close(10)

        !输出最优汇流参数
        filename=trim(address)//'\结果\最优汇流参数.txt'
        open(10,file=filename)
        do i=1,NRR
            write(10,*)Xbest(NRP+i)
        end do
        close(10)

        !输出每场洪水、每个雨量站对应的最优状态量
        if(PRModel==1)then
            N1=5
        else
            N1=2
        end if
        filename=trim(address)//'\结果\最优状态量.txt'
        open(10,file=filename)        
        do i=1,NF
            write(Number,"(I4)")i
            Cline='第'//trim(adjustl(Number))//'场洪水'
            write(10,*)trim(adjustl(Cline))
            do j=1,NP
                write(10,*)Ybest((i-1)*NP*N1+(j-1)*N1+1:(i-1)*NP*N1+j*N1)
            end do
        end do
        close(10)

        !输出每场洪水的模拟结果
        do i=1,NF
            write(Number,"(I4)")i
            filename=trim(Address)//'\结果\第'//trim(adjustl(Number))//'场洪水模结果.txt'
            open(10,file=filename)
            write(10,*)'实测值 模拟值'
            do j=1,FloodTime(i+1)-FloodTime(i)
               write(10,*)Q(FloodTime(i)+j),QC(FloodTime(i)+j)
            end do
            close(10)
        end do

    end subroutine
end module