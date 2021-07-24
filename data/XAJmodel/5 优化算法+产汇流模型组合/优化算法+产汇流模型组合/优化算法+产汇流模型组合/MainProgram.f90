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
    filename=trim(address)//'\�������ļ�\��ʼ����.txt'
    open(10,file=filename)
    read(10,*)AA                       !�������
    read(10,*)PRModel                  !����ģ��
    read(10,*)RoutingModel             !����ģ��
    read(10,*)CarModel                 !�����ʶ��㷨
    read(10,*)NF                       !��ˮ����
    read(10,*)NP                       !����վ����
    read(10,*)Key1                     !Key1Ϊ1ʱ״̬���Ѹ�����2ʱ���ʶ�״̬��
    read(10,*)Key2                     !Key2Ϊ1ʱ��ģ�Ͳ����ʶ���2ʱΪԤ��
    read(10,*)Iter                     !�Ż��㷨ѭ������
    read(10,*)Random_On                !�Ƿ���������(����ÿ���Ż�����᲻ͬ��������ÿ��һ��)��1ʱΪ��������ԣ�2ʱΪȥ�������
    read(10,*)IN                       !���ο���վ����
    close(10)
    
    call readdata(1)                   !��������
    call ParameterScale(Key2)          !������Χ����ָ����


    !�����Ż���CarModelΪ1ʱ����SCE-UA�㷨��Ϊ2ʱ���õ����Է�
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

    !������
    call DataOutput(QC)
end program