module Calibration_SM
    use global
    use ModelCombination
    implicit none

    contains
    subroutine SM(upbound,downbound,xx,OFV,QC)
      !µ¥´¿ÐÔ·¨

       implicit none
       real::upbound(:),downbound(:),xx(:)
       real,allocatable::QC(:),ff(:)
       real::OFV
       integer::n,cato,nf,periods,sz,nn
       integer::t=100000
       integer::k,i,j,iterations,status,maxfloc,minfloc
       real::b,a,d1,d2,f,f_reflect,f_exte,f_xc,maxf,minf
       real::x(size(upbound)+1,size(upbound)+1),center(size(upbound)),reflect(size(upbound)),exte(size(upbound)),xc(size(upbound))
       real,allocatable::r(:)
       integer::x_up(size(upbound)),x_down(size(upbound)),y(1)

       nn=0
       b=0.1;a=1.0
       n=size(upbound)
       If(Random_On==1)then
           call random_seed()
       end if
       allocate(r(n))
       call random_number(r(1:n))
       x(1,1:n)=r*(upbound-downbound)+downbound
       allocate(ff(10000))
       if(Key3==1)then
           call CombinedModel(x(1,1:n),Ybest,QC,x(1,n+1))
       else
           call CombinedModel(Xbest,x(1,1:n),QC,x(1,n+1))
       end if
       OFV=x(1,n+1)
       d1=((n+1)**0.5+n-1)*b/2.0**0.5/n
       d2=((n+1)**0.5-1)*b/2.0**0.5/n

       do i=2,n+1
          x(i,1:n)=x(1,1:n)
          x(i,1:n)=x(i,1:n)+d2
          x(i,i-1)=x(i,i-1)-d2+d1
          if(Key3==1)then
              call CombinedModel(x(i,1:n),Ybest,QC,x(i,n+1))
          else
              call CombinedModel(Xbest,x(i,1:n),QC,x(i,n+1))
          end if           
       end do

       do iterations=1,t
          y=maxloc(x(:,n+1))
          maxfloc=y(1)
          maxf=maxval(x(:,n+1))
          y=minloc(x(:,n+1))
          minfloc=y(1)
          minf=minval(x(:,n+1))

          center=(sum(x(:,1:n),1)-x(maxfloc,1:n))/n
          reflect=(1+a)*center-a*x(maxfloc,1:n)

          x_down(1:n)=(/(0,i=1,n)/)
          x_up(1:n)=(/(0,i=1,n)/)

          do j=1,n
            if (upbound(j)<reflect(j))   x_up(j)=1  
            if (downbound(j)>reflect(j))  x_down(j)=1
          end do

          call random_seed()
          call random_number(r)
          reflect=(reflect-reflect*x_up)+x_up*center+r*(x_up*upbound-x_up*center)

          call random_seed()
          call random_number(r)
          reflect=(reflect-reflect*x_down)+x_down*center-r*(x_down*center-x_down*downbound)

          if(Key3==1)then
              call CombinedModel(Reflect,Ybest,QC,f_reflect)
          else
              call CombinedModel(Xbest,Reflect,QC,f_reflect)
          end if 
          
          if (f_reflect<maxf) then
              exte=(1+2.0*a)*center-2.0*a*x(maxfloc,1:n)

              x_down(1:n)=(/(0,i=1,n)/)
              x_up(1:n)=(/(0,i=1,n)/)

              do j=1,n
                if (upbound(j)<exte(j))   x_up(j)=1 
                if (downbound(j)>exte(j))  x_down(j)=1
              end do
              
              call random_seed()
              call random_number(r)
              exte=(exte-exte*x_up)+x_up*center+r*(x_up*upbound-x_up*center)

              call random_seed()
              call random_number(r)
              exte=(exte-exte*x_down)+x_down*center-r*(x_down*center-x_down*downbound)

              if(Key3==1)then
                  call CombinedModel(exte,Ybest,QC,f_exte)
              else
                  call CombinedModel(Xbest,exte,QC,f_exte)
              end if        
  
              if (f_exte<maxf) then
                  x(maxfloc,1:n)=exte
                  x(maxfloc,n+1)=f_exte
              else
                  x(maxfloc,1:n)=reflect
                  x(maxfloc,n+1)=f_reflect
              end if

         else
             xc=center-(center-x(maxfloc,1:n))*a/2.0
             
             x_down(1:n)=(/(0,i=1,n)/)
             x_up(1:n)=(/(0,i=1,n)/)

             do j=1,n
               if (upbound(j)<xc(j))   x_up(j)=1  
               if (downbound(j)>xc(j))  x_down(j)=1
             end do

             call random_seed()
             call random_number(r)
             xc=(xc-xc*x_up)+x_up*center+r*(x_up*upbound-x_up*center)

             call random_seed()
             call random_number(r)
             xc=(xc-xc*x_down)+x_down*center-r*(x_down*center-x_down*downbound)

             if(Key3==1)then
                 call CombinedModel(xc,Ybest,QC,f_xc)
             else
                 call CombinedModel(Xbest,xc,QC,f_xc)
             end if 

             if (f_xc<maxf) then
                x(maxfloc,1:n)=xc
                x(maxfloc,n+1)=f_xc
             else
                do k=1,n+1
                   if (k/=minfloc) then
                      x(k,1:n)=(x(minfloc,1:n)+x(k,1:n))/2.0

                      x_down(1:n)=(/(0,i=1,n)/)
                      x_up(1:n)=(/(0,i=1,n)/)
                      do j=1,n
                         if (upbound(j)<x(k,j))   x_up(j)=1  
                         if (downbound(j)>x(k,j))  x_down(j)=1
                      end do
                      
                      call random_seed()
                      call random_number(r)
                      x(k,1:n)=(x(k,1:n)-x(k,1:n)*x_up)+x_up*x(k,1:n)+r*(x_up*upbound-x_up*x(k,1:n))

                      call random_seed()
                      call random_number(r)
                      x(k,1:n)=(x(k,1:n)-x(k,1:n)*x_down)+x_down*x(k,1:n)-r*(x_down*x(k,1:n)-x_down*downbound)

                     if(Key3==1)then
                         call CombinedModel(x(k,1:n),Ybest,QC,x(k,n+1))
                     else
                         call CombinedModel(Xbest,x(k,1:n),QC,x(k,n+1))
                     end if
                  end if

                end do
            end if
         end if 
         if(abs(mod(iterations,20))<1E-4)then
           nn=nn+1
           y=minloc(x(:,n+1))
           minfloc=y(1)
           maxf=x(minfloc,n+1) 
           write(*,*)'iterations=',iterations
           write(*,*)maxf
           ff(nn)=maxf
         end if    
         if(iterations>400)then
           if(abs(sum(ff(nn-10:nn-1))/10-ff(nn))<1E-5)then
             exit
           end if
         end if
       end do 
       y=minloc(x(:,n+1))
       minfloc=y(1);xx=x(minfloc,1:n)
 
    end subroutine    


end module