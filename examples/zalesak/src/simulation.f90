!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Only get a VOF solver and corresponding time tracker
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Velocity arrays
   real(WP), dimension(:,:,:), allocatable :: U,V,W
   real(WP), dimension(3) :: center
   real(WP) :: radius,width,height
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Initialize our VOF
      initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: VFhi,VFlo,remap,lvira,plicnet,r2pnet,r2p
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=r2pnet,transport_method=remap,name='VOF')
         ! Initialize to Zalesak disk
         center=[0.0_WP,0.25_WP,0.0_WP]
         radius=0.05_WP!0.15_WP
         width =0.05_WP
         height=0.25_WP
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  !call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_zalesak,0.0_WP,amr_ref_lvl)
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volume
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block initialize_vof
      
      
      ! Prepare our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         integer :: i,j,k
         ! Allocate arrays
         allocate(U(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_))
         allocate(V(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_))
         allocate(W(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_))
         ! Initialize to analytical flow field
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Solid body rotation about z-axis
                  !U(i,j,k)=-twoPi*vf%cfg%ym(j)
                  !V(i,j,k)=+twoPi*vf%cfg%xm(i)
                  !W(i,j,k)=0.0_WP
                  ! Radially symmetric stagnation flow
                  U(i,j,k)=+1.0_WP*vf%cfg%x(i)
                  V(i,j,k)=-2.0_WP*vf%cfg%y(j)
                  W(i,j,k)=+1.0_WP*vf%cfg%z(k)
               end do
            end do
         end do
      end block initialize_velocity
      
      
      ! Initialize time tracker - only 1 sub-iteration!
      initialize_timetracker: block
         time=timetracker(amRoot=vf%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         integer :: i,j,k,np,nplane
         smesh=surfmesh(nvar=2,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='thickness'
         ! Transfer polygons to smesh
         call vf%update_surfmesh_nowall(smesh)
         ! Calculate thickness
         call vf%get_thickness()
         ! Populate nplane and thickness variables
         smesh%var(1,:)=1.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                        smesh%var(2,np)=vf%thickness(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='zalesak')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call vf%get_max()
         call vf%get_cfl(dt=time%dt,U=U,V=V,W=W,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(vf%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%write()
      end block create_monitor
      
      
   contains
      
      
      !> Function that defines a level set function for Zalesak's notched circle problem
      function levelset_zalesak(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         real(WP) :: c,b,b1,b2,h1,h2
         c=radius-sqrt(sum((xyz-center)**2))
         b1=center(1)-0.5_WP*width; b2=center(1)+0.5_WP*width
         h1=center(2)-radius*cos(asin(0.5_WP*width/radius)); h2=center(2)-radius+height
         if     (c.ge.0.0_WP.and.xyz(1).le.b1.and.xyz(2).le.h2) then
            b=b1-xyz(1); G=min(c,b)
         else if (c.ge.0.0_WP.and.xyz(1).ge.b2.and.xyz(2).le.h2) then
            b=xyz(1)-b2; G=min(c,b)
         else if (c.ge.0.0_WP.and.xyz(1).ge.b1.and.xyz(1).le.b2.and.xyz(2).ge.h2) then
            b=xyz(2)-h2; G=min(c,b)
         else if (c.ge.0.0_WP.and.xyz(1).le.b1.and.xyz(2).ge.h2) then
            b=sqrt(sum((xyz-(/b1,h2,0.0_WP/))**2)); G=min(c,b)
         else if (c.ge.0.0_WP.and.xyz(1).ge.b2.and.xyz(2).ge.h2) then
            b=sqrt(sum((xyz-(/b2,h2,0.0_WP/))**2)); G=min(c,b)
         else if (xyz(1).ge.b1.and.xyz(1).le.b2.and.xyz(2).le.h2.and.xyz(2).ge.h1) then
            G=-min(abs(xyz(1)-b1),abs(xyz(1)-b2),abs(xyz(2)-h2))
         else if (xyz(1).ge.b1.and.xyz(1).le.b2.and.xyz(2).le.h1) then
            G=-min(sqrt(sum((xyz-(/b1,h1,0.0_WP/))**2)),sqrt(sum((xyz-(/b2,h1,0.0_WP/))**2)))
         else
            G=c
         end if
      end function levelset_zalesak
      
      
      !> Function that defines a level set function for a sphere
      function levelset_sphere(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=radius-sqrt(sum((xyz)**2))
      end function levelset_sphere
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call vf%get_cfl(dt=time%dt,U=U,V=V,W=W,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ================ VOF SOLVER =======================
            call vf%advance(dt=time%dt,U=U,V=V,W=W)
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_smesh: block
               use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
               integer :: i,j,k,np,nplane
               ! Transfer polygons to smesh
               call vf%update_surfmesh_nowall(smesh)
               ! Calculate thickness
               call vf%get_thickness()
               ! Also populate nplane variable
               smesh%var(1,:)=1.0_WP
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                              smesh%var(2,np)=vf%thickness(i,j,k)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call vf%get_max()
         call mfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(U,V,W)
   end subroutine simulation_final
   
   
end module simulation