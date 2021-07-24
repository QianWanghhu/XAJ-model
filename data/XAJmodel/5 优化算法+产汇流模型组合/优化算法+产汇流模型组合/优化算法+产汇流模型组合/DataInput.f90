module DataInput
    use global
    implicit none

    contains
    subroutine readdata(n)
        implicit none
        integer::i,j,IS,n,N1

        !���뽵ˮ����������������
        filename=trim(Address)//'\�������ļ�\data.txt'
        open(10,file=filename)
        TT=0
        do
		    read(10,*,IOSTAT=IS)
		    if(IS/=0) exit
		    TT=TT+1
	   end do
       TT=TT-1
       allocate(P(TT,NP))
       allocate(Q(TT))
       allocate(E(TT))
       rewind(10)
       read(10,*)
       do i=1,TT
           read(10,*)E(i),Q(i),P(i,:)
       end do
       close(10)

        !���뽵ˮ����������������
        if(IN>0)then
            filename=trim(Address)//'\�������ļ�\���ο���վ.txt'
            open(10,file=filename)
            allocate(QIN(TT,IN))
            read(10,*)
            do i=1,TT
                read(10,*)QIN(i,:)
            end do
            close(10)
        end if


       !����վȨ��
       allocate(weight(NP)) 
       filename=trim(Address)//'\�������ļ�\����վȨ��.txt'
       open(10,file=filename)
       do i=1,NP
            read(10,*)weight(i)
       end do
       close(10) 

       allocate(FloodTime(NF+1))
       FloodTime(1)=0
       filename=trim(Address)//'\�������ļ�\�κ�ʱ��.txt'
       open(10,file=filename)
       read(10,*)
       do i=1,NF
           read(10,*)FloodTime(i+1)
       end do
       close(10)

       if(PRModel==1)then
           N1=5
       else
           N1=2
       end if
       allocate(state(NF,NP,N1),Ybest(NF*NP*N1))
       state=0.0;Ybest=0.01

       !��ʼ״̬��
       if(Key1==1)then
           filename=trim(Address)//'\�������ļ�\��ʼ״̬��.txt'
           open(10,file=filename)
           read(10,*)
           do i=1,NF
               do j=1,NP
                   read(10,*)state(i,j,:)
                   Ybest((i-1)*NP*N1+(j-1)*N1+1:(i-1)*NP*N1+j*N1)=state(i,j,:)
               end do
           end do
           close(10) 
       end if

       !Ŀ�꺯��Ȩ��
       filename=trim(Address)//'\�������ļ�\Ŀ�꺯������.txt'
       open(10,file=filename)
       do i=1,3
           read(10,*)CC(i)
       end do
       close(10)
    end subroutine
end module