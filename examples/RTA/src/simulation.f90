!> Various definitions and tools for running an NGA2 simulation
module simulation
   use hit_class, only: hit
   use rta_class, only: rta
   implicit none
   private
   
   !> HIT simulation
   type(hit) :: turb
   logical :: isInHITGrp

   !> Rayleigh-Taylor atomization simulation
   type(rta) :: atom
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: hit_group
      
      ! Initialize atomization simulation
      call atom%init()
      
      ! Create an MPI group using leftmost processors only
      create_hit_group: block
         use parallel, only: group,comm
         use mpi_f08,  only: MPI_Group_incl
         integer, dimension(:), allocatable :: ranks
         integer, dimension(3) :: coord
         integer :: n,ngrp,ierr,ny,nz
         ngrp=atom%cfg%npy*atom%cfg%npz
         allocate(ranks(ngrp))
         ngrp=0
         do nz=1,atom%cfg%npz
            do ny=1,atom%cfg%npy
               ngrp=ngrp+1
               coord=[0,ny-1,nz-1]
               call MPI_CART_RANK(atom%cfg%comm,coord,ranks(ngrp),ierr)
            end do
         end do
         call MPI_Group_incl(group,ngrp,ranks,hit_group,ierr)
         if (atom%cfg%iproc.eq.1) then
            isInHITGrp=.true.
         else
            isInHITGrp=.false.
         end if
      end block create_hit_group
      
      ! Initialize HIT simulation
      if (isInHITGrp) call turb%init(group=hit_group)
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      do while (.not.atom%time%done())
         
         ! Advance atomization simulation
         call atom%step()
         
         ! Advance HIT simulation and transfer velocity info
         if (isInHITGrp) then
            ! Advance HIT
            call turb%step(atom%time%dt)
            ! Transfer turbulent velocity from hit to rta
            apply_boundary_condition: block
               use tpns_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,ihit
               call atom%fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  ihit=i-atom%fs%cfg%imin+turb%fs%cfg%imax+1
                  atom%fs%U(i  ,j,k)=turb%fs%U(ihit  ,j,k)+atom%Uconv
                  atom%fs%V(i-1,j,k)=turb%fs%V(ihit-1,j,k)
                  atom%fs%W(i-1,j,k)=turb%fs%W(ihit-1,j,k)
               end do
            end block apply_boundary_condition
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize atomization simulation
      call atom%final()
      
      ! Finalize HIT simulation
      if (isInHITGrp) call turb%final()
      
   end subroutine simulation_final
   

end module simulation