module Data_Output
    use global
    implicit none

    contains
    !������
    subroutine DataOutput(QC)
        implicit none
        character(len=2000)::Number,Cline
        real::QC(:)
        integer::i,j,N1
        
        !������Ų�������
        filename=trim(address)//'\���\���Ų�������.txt'
        open(10,file=filename)
        do i=1,NRP
            write(10,*)Xbest(i)
        end do
        close(10)

        !������Ż�������
        filename=trim(address)//'\���\���Ż�������.txt'
        open(10,file=filename)
        do i=1,NRR
            write(10,*)Xbest(NRP+i)
        end do
        close(10)

        !���ÿ����ˮ��ÿ������վ��Ӧ������״̬��
        if(PRModel==1)then
            N1=5
        else
            N1=2
        end if
        filename=trim(address)//'\���\����״̬��.txt'
        open(10,file=filename)        
        do i=1,NF
            write(Number,"(I4)")i
            Cline='��'//trim(adjustl(Number))//'����ˮ'
            write(10,*)trim(adjustl(Cline))
            do j=1,NP
                write(10,*)Ybest((i-1)*NP*N1+(j-1)*N1+1:(i-1)*NP*N1+j*N1)
            end do
        end do
        close(10)

        !���ÿ����ˮ��ģ����
        do i=1,NF
            write(Number,"(I4)")i
            filename=trim(Address)//'\���\��'//trim(adjustl(Number))//'����ˮģ���.txt'
            open(10,file=filename)
            write(10,*)'ʵ��ֵ ģ��ֵ'
            do j=1,FloodTime(i+1)-FloodTime(i)
               write(10,*)Q(FloodTime(i)+j),QC(FloodTime(i)+j)
            end do
            close(10)
        end do

    end subroutine
end module