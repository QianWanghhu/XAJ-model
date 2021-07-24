!全局变量模块
module XAJ_global
    implicit none  
    real,allocatable,save::Q(:),E(:),P(:)            
    character(len=2000),save::Address,filename
end module