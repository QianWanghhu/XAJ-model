program main
    use XAJ_global
    use XAJ_RunoffProduction
    implicit none

    real,allocatable::x(:),Output(:,:),y(:)  
    integer::i,IS,szz,j
    real::mu,AA,AA1

    call getcwd(Address)

    !读入降水、蒸发、径流数据
    filename=trim(Address)//'\需设置文件\data.txt'
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

   !读入模型参数
   allocate(x(10)) 
   filename=trim(Address)//'\需设置文件\模型参数.txt'
   open(10,file=filename)
   do i=1,10
      read(10,*)x(i)
   end do
   close(10)

   !读入初始状态量
   allocate(y(5))
   filename=trim(Address)//'\需设置文件\初始状态量.txt'
   open(10,file=filename)
   do i=1,5
       read(10,*)y(i)
   end do
   close(10)  
     
   !计算每个雨量站所在区域的产流过程
   allocate(Output(szz,3))
   Output=0.0
   Call XAJ_RP(x,y,24,Output)

   !输出整个流域的产流过程,分为三列输出，从左到右分别为地表径流、壤中流和地下径流
   filename=trim(Address)//'\结果\产流总过程.txt'
   open(20,file=filename)
   write(20,*)'地表产流 壤中流产流量 地下产流量'
   do i=1,szz
      write(20,"(3f10.5)")Output(i,:)
   end do
   close(20)

end program