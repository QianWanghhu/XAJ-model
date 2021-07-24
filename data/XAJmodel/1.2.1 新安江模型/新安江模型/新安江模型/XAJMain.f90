program main
    use XAJ_global
    use XAJ_RunoffProduction
    implicit none

    real,allocatable::x(:),Output(:,:),y(:)  
    integer::i,IS,szz,j
    real::mu,AA,AA1

    call getcwd(Address)

    !���뽵ˮ����������������
    filename=trim(Address)//'\�������ļ�\data.txt'
    open(10,file=filename)
    szz=0
	do
		read(10,*,IOSTAT=IS)
		if(IS/=0) exit
		szz=szz+1
    end do
    szz=szz-1
    allocate(P(szz))
    allocate(Q(szz))
    allocate(E(szz))
    rewind(10)
    read(10,*)
    do i=1,szz
       read(10,*)P(i),E(i),Q(i)
    end do
    close(10)

   !����ģ�Ͳ���
   allocate(x(10)) 
   filename=trim(Address)//'\�������ļ�\ģ�Ͳ���.txt'
   open(10,file=filename)
   do i=1,10
      read(10,*)x(i)
   end do
   close(10)

   !�����ʼ״̬��
   allocate(y(5))
   filename=trim(Address)//'\�������ļ�\��ʼ״̬��.txt'
   open(10,file=filename)
   do i=1,5
       read(10,*)y(i)
   end do
   close(10)  
     
   !����ÿ������վ��������Ĳ�������
   allocate(Output(szz,3))
   Output=0.0
   Call XAJ_RP(x,y,24,Output)

   !�����������Ĳ�������,��Ϊ��������������ҷֱ�Ϊ�ر������������͵��¾���
   filename=trim(Address)//'\���\�����ܹ���.txt'
   open(20,file=filename)
   write(20,*)'�ر���� ������������ ���²�����'
   do i=1,szz
      write(20,"(3f10.5)")Output(i,:)
   end do
   close(20)

end program