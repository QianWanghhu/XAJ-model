program Carlibration_Main
    use global
    use Carlibration_SCE
    use Calibration_SM
    use DataInput
    use Data_Output

    implicit none
    character(len=2000)::Number,Cline
    integer::i,j,Iter
    real::OFV
    real,allocatable::QC(:)

    call getcwd(address)
    filename=trim(address)//'\需设置文件\初始条件.txt'
    open(10,file=filename)
    read(10,*)AA                       !流域面积
    read(10,*)PRModel                  !产流模型
    read(10,*)RoutingModel             !汇流模型
    read(10,*)CarModel                 !参数率定算法
    read(10,*)NF                       !洪水场次
    read(10,*)NP                       !雨量站个数
    read(10,*)Key1                     !Key1为1时状态量已给出，2时需率定状态量
    read(10,*)Key2                     !Key2为1时是模型参数率定，2时为预测
    read(10,*)Iter                     !优化算法循环次数
    read(10,*)Random_On                !是否加入随机性(加入每次优化结果会不同，不加入每次一样)：1时为加入随机性，2时为去除随机性
    read(10,*)IN                       !上游控制站个数
    close(10)
    
    call readdata(1)                   !数据输入
    call ParameterScale(Key2)          !参数范围或是指输入


    !参数优化，CarModel为1时采用SCE-UA算法，为2时采用单纯性法
    allocate(QC(TT))
    do i=1,iter
        if(CarModel==1)then        
            if(Key2==1)then
                Key3=1
                call SCE_UA(UpB,DownB,Xbest,OFV,QC)
                call StateBoundGet(Xbest,UpState,DownState)
                if(Key1==2)Then
                    key3=2
                    call SCE_UA(UpState,DownState,Ybest,OFV,QC)
                End If
            else
                key3=2
                call StateBoundGet(Xbest,UpState,DownState)
                call SCE_UA(UpState,DownState,Ybest,OFV,QC)
            end if
        else
            if(Key2==1)then
                Key3=1
                call SM(UpB,DownB,Xbest,OFV,QC)
                call StateBoundGet(Xbest,UpState,DownState)
                if(Key1==2)then
                    Key3=2
                    call SM(UpState,DownState,Ybest,OFV,QC)
                end if
            else
                key3=2
                call StateBoundGet(Xbest,UpState,DownState)
                call SM(UpState,DownState,Ybest,OFV,QC)
            end if
        end if
    end do

    !结果输出
    call DataOutput(QC)
end program