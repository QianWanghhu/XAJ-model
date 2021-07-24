module global

    implicit none
    real,allocatable,save::Q(:),E(:),P(:,:),QIN(:,:)
    real,allocatable,save::weight(:),state(:,:,:)
    integer::TT,NP,NF,Key1,Key2,PRModel,RoutingModel,CarModel,Key3,Random_On,IN
    integer,allocatable::FloodTime(:)
    Character(len=2000)::filename,address
    real::AA
    real,allocatable::UpB(:),DownB(:),UpState(:),DownState(:),Xbest(:),Ybest(:)
    integer::NRP,NRR
    real::periods=1.0
    real::CC(3)

end module