!> Various definitions and tools for running an NGA2 simulation
module simulation
   use counterflow_class, only: counterflow
   implicit none
   private
   
   !> Counterflow jet simulation
   type(counterflow) :: cflow
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize counterflow simulation
      call cflow%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Counterflow drives overall time integration
      do while (.not.cflow%time%done())
         
         ! Advance counterflow simulation
         call cflow%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize counterflow simulation
      call cflow%final()
      
   end subroutine simulation_final
   
   
end module simulation