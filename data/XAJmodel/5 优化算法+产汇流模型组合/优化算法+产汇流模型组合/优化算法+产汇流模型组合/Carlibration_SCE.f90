module Carlibration_SCE
    use global
    use ModelCombination
    implicit none

    contains
    
    !SCE-UA�㷨
    subroutine SCE_UA(upbound,downbound,xx,OFV,QC)
        implicit none
        !���ó�ʼ����
        integer::n,m,s,qq,o,alp,beta,dd,pp
        integer::i,j,k,kk,SCE,cato,model
        real::rnd2,maxf,OFV
        integer::sample,L0,t,sz,nf
        real::upbound(:),downbound(:),xx(:),QC(:)
        real,allocatable::rnd(:)
        real,allocatable::x1(:,:)
        real,allocatable::f(:)
        real,allocatable::D(:,:)
        real,allocatable::A(:,:,:)
        real,allocatable::C(:,:)
        real,allocatable::L(:)
        integer,allocatable::Lx(:)
        integer,allocatable::counts(:)
        integer,allocatable::Cou(:)
        real,allocatable::b(:,:)
        real,allocatable::B0(:,:)
        integer,allocatable::BF(:)
        real,allocatable::center(:)
        real,allocatable::Reflect(:)
        real,allocatable::Up(:)
        real,allocatable::Low(:)
        real::f_Ref,f_Contra,min_n,threshold,mu,f_Mua
        integer::min_l
        real,allocatable::Contra(:)
        real,allocatable::Mutation(:)
        
        n=size(upbound)             !���Ż��Ĳ�������
        m=n+1                       !�������е�ĸ���
        pp=2                         !�����εĸ���
        s=m*pp                       !������
        qq=3                       !ÿ���ӵ������е�ĸ���
        o=100                        !SCE-UA�㷨��ִ�д���
        alp=1                       !�����������Ӵ�����
        beta=2*n+1                  !ÿ�������ν�������Ĳ�����
        dd=1
        allocate(rnd(n))    
        allocate(x1(s,n+1))
        allocate(D(s,n+1))
        allocate(f(s))
        allocate(counts(m))
        allocate(A(m,n+1,pp))
        allocate(C(m,n+1))
        allocate(L(qq))
        allocate(Cou(m))
        allocate(b(qq,n+1))
        allocate(B0(qq,n+1))
        allocate(BF(qq))
        allocate(center(n))
        allocate(Reflect(n))
        allocate(Up(n))
        allocate(Low(n))
        allocate(Contra(n))
        allocate(Mutation(n))


        !�ں���Ķ������������ȡs�������
        If(Random_On==1)then
            call random_seed()
        end if
        do i=1,s
            call random_number(rnd)
            x1(i,1:n)=rnd*(upbound-downbound)+downbound
            if(Key3==1)then
                call CombinedModel(x1(i,1:n),Ybest,QC,OFV)
            else
                call CombinedModel(Xbest,x1(i,1:n),QC,OFV)
            end if
            f(i)=OFV
            x1(i,n+1)=f(i)
        end do
                
        !����Ŀ�꺯���Բ�����Ӵ�С����
        D=x1
        allocate(Lx(size(D,1)))
        call sortrows(D,n+1,Lx)
        counts(1)=m
        do i=2,m
            counts(i)=counts(i-1)+m-i+1
        end do

        !SCE-UA�㷨����
        do SCE=1,o
            !����Ŀ�꺯���Բ������С��������
            call sortrows(D,n+1,Lx)
            !��s�����������pp���Ӿ���
            do i=1,pp
                do j=1,m
                    A(j,:,i)=D(i+pp*(j-1),:)
                end do
            end do
            do i=1,Beta
                do j=1,pp
                    C=A(:,:,j)
                    L=0
                    Cou=counts
                    !���ݸ��ʣ����Ӿ����г�ȡqq��������
                    do k=1,qq
                        call random_number(rnd2)
                        sample=ceiling(rnd2*maxval(Cou))
                        L0=1
                        do kk=1,m
                            if (Cou(kk)>=sample) then
                                L0=kk
                                exit
                            end if
                        end do
                        L(k)=L0
                        do kk=1,m
                            if (Cou(kk)>=sample) then
                                Cou(kk)=Cou(kk)-(m-L(k)+1)
                            end if
                        end do
                        b(k,:)=C(L(k),:)                   !��ѡȡ��q�������b��
                    end do
                    !�������Ż�
                    do k=1,alp
                        B0=b                               !����sortrows��b��ı�
                        call sortrows(b,n+1,BF)
                        center=(sum(B0(:,1:n),1)-B0(BF(qq),1:n))/(qq-1)   !�����ĵ�
                        Reflect=center+dd*(center-B0(BF(qq),1:n))         !�㷴���
                        Up=maxval(B0(:,1:n),1)
                        Low=minval(B0(:,1:n),1)
                        if ((.not. all(upbound>Reflect)) .or. (.not. all(Reflect>downbound))) then
                            call random_number(rnd)
                            Reflect=rnd*(Up-Low)+Low
                        end if 

                        if(Key3==1)then
                            call CombinedModel(Reflect,Ybest,QC,OFV)
                        else
                            call CombinedModel(Xbest,Reflect,QC,OFV)
                        end if
                                                                    
                        f_Ref=OFV
                        if (f_Ref<B0(BF(qq),n+1)) then
                            B0(BF(qq),1:n)=Reflect
                            B0(BF(qq),n+1)=f_Ref
                        else
                            Contra=center-(center-B0(BF(qq),1:n))*dd/2

                            if(Key3==1)then
                                call CombinedModel(Contra,Ybest,QC,OFV)
                            else
                                call CombinedModel(Xbest,Contra,QC,OFV)
                            end if

                            f_Contra=OFV
                            if (f_Contra<B0(BF(qq),n+1)) then
                                B0(BF(qq),1:n)=Contra
                                B0(BF(qq),n+1)=f_Contra
                            else
                                call random_number(rnd)
                                Mutation=rnd*(Up-Low)+Low

                                if(Key3==1)then
                                    call CombinedModel(Mutation,Ybest,QC,OFV)
                                else
                                    call CombinedModel(Xbest,Mutation,QC,OFV)
                                end if     

                                f_Mua=OFV
                                B0(BF(qq),1:n)=Mutation
                                B0(BF(qq),n+1)=f_Mua
                            end if
                        end if
                        b=B0
                    end do

                    !���ص��Ӿ���
                    do k=1,qq
                        C(L(k),:)=b(k,:)
                    end do 

                    !����
                    call sortrows(C,n+1,Lx)
                    A(:,:,j)=C
                end do
            end do

            !�������η����Ż��������ԭ����D��
            do i=1,pp
                do j=1,m
                    D(i+pp*(j-1),:)=A(j,:,i)
                end do
            end do
            min_n=minval(D(:,n+1))
            min_l=0
            do i=1,n+1
                min_l=min_l+1
                if(min_n==D(i,n+1)) exit
            end do
            threshold=1-min_n
            write(*,*)'SCE=',SCE
            write(*,*)'OFV=',min_n
        end do
        xx=D(min_l,1:n);maxf=D(min_l,n+1)      
               
    end subroutine
!--------------------------------------------------------------------------------------------!
    !�Ծ������������򣬴�С����
    !D��������ľ���nn��������յ��У�Lx����������ԭ�к�
    !r0:����������c0����������
    subroutine sortrows(D,nn,Lx)
        implicit none
        real::D(:,:)
        integer::Lx(:)
        integer::nn,r0,c0
        integer::i,j
        real,allocatable::temp(:)
        real,allocatable::D_1(:,:)          
        r0=size(D,1)
        c0=size(D,2)
        allocate(temp(c0))
        allocate(D_1(r0,c0))
        D_1=D                               !�洢ԭʼ������D
        do i=1,r0-1
            do j=i+1,r0
                if(D(i,nn)>D(j,nn))then
                    temp=D(i,:)
                    D(i,:)=D(j,:)
                    D(j,:)=temp
                end if
            end do
        end do
        do i=1,r0
            do j=1,r0
                if (all(abs(D(i,:)-D_1(j,:))<1.0e-34))then
                    Lx(i)=j
                    exit
                end if
            end do
        end do
    end subroutine
end module