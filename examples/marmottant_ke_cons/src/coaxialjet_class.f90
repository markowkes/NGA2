!> Definition for a coaxial jet class
module coaxialjet_class
   use precision,         only: WP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   use stracker_class,    only: stracker
   implicit none
   private
   
   public :: coaxialjet
   
   !> roundcoaxialjet jet object
   type :: coaxialjet
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf     !< Volume fraction solver
      type(tpns)        :: fs     !< Two-phase flow solver
      type(hypre_str)   :: ps     !< Structured Hypre linear solver for pressure
      type(timetracker) :: time   !< Time info
      
      !> Implicit velocity solver
      type(ddadi) :: vs           !< DDADI solver for velocity
      
      !> SGS modeling
      logical        :: use_sgs   !< Is an LES model used?
      type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh     !< Surface mesh for interface
      type(ensight)  :: ens_out   !< Ensight output for flow variables
      type(event)    :: ens_evt   !< Event trigger for Ensight output
      
      !> Simulation monitoring files
      type(monitor) :: mfile      !< General simulation monitoring
      type(monitor) :: cflfile    !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      integer  :: nlayer=4                 !< Size of buffer layer for VOF removal
      
      !> Provide a pardata object for restarts
      logical       :: restarted  !< Is the simulation restarted?
      type(pardata) :: df         !< Pardata object for restart I/O
      type(event)   :: save_evt   !< Event to trigger restart I/O

      !> Flow definition
      real(WP) :: Dl,Hg,Ul,Ug,dg,lip
      
   contains
      
      procedure :: init                            !< Initialize round jet simulation
      procedure :: step                            !< Advance round jet simulation by one time step
      procedure :: final                           !< Finalize round jet simulation
      
   end type coaxialjet
   
   
contains
   
   
   !> Initialization of roundjet simulation
   subroutine init(this)
      use param, only: param_read
      implicit none
      class(coaxialjet), intent(inout) :: this
      
      
      ! Initialize the config
      initialize_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x_uni(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y_uni(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z_uni(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! Add stretching
         call param_read('Stretched cells in yz',ns_yz,default=0)
         if (ns_yz.gt.0) call param_read('Stretch ratio in yz',sratio_yz)
         call param_read('Stretched cells in x' ,ns_x ,default=0)
         if (ns_x .gt.0) call param_read('Stretch ratio in x' ,sratio_x )
         allocate(x(nx+1+1*ns_x )); x(      1:      1+nx)=x_uni
         allocate(y(ny+1+2*ns_yz)); y(ns_yz+1:ns_yz+1+ny)=y_uni
         allocate(z(nz+1+2*ns_yz)); z(ns_yz+1:ns_yz+1+nz)=z_uni
         do i=nx+2,nx+1+ns_x
            x(i)=x(i-1)+sratio_x*(x(i-1)-x(i-2))
         end do
         do j=ns_yz,1,-1
            y(j)=y(j+1)+sratio_yz*(y(j+1)-y(j+2))
         end do
         do j=ns_yz+2+ny,ny+1+2*ns_yz
            y(j)=y(j-1)+sratio_yz*(y(j-1)-y(j-2))
         end do
         do k=ns_yz,1,-1
            z(k)=z(k+1)+sratio_yz*(z(k+1)-z(k+2))
         end do
         do k=ns_yz+2+nz,nz+1+2*ns_yz
            z(k)=z(k-1)+sratio_yz*(z(k-1)-z(k-2))
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='jet')
         ! Read in partition
         call param_read('Partition',partition)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         ! No walls in the atomization domain
         this%cfg%VF=1.0_WP
      end block initialize_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resU         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,plicnet,flux
         use mms_geom,  only: cube_refine_vol
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with plicnet reconstruction
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')
         ! Initialize to clipped cylindrical interface
         call param_read('Liquid diameter',this%Dl)
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cylinder,0.0_WP,amr_ref_lvl)
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  else
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
                  ! Clip cylinder after two cells
                  if (i.ge.this%cfg%imin+2) then
                     this%vf%VF(i,j,k)=0.0_WP
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set simple full-liquid/full-gas interface planes in geometric overlap cells
         call this%vf%set_full_bcond()
         ! Now apply Neumann condition on interface at inlet to have proper round injection
         neumann_irl: block
            use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
            &                                setNumberOfPlanes,setPlane,matchVolumeFraction
            real(WP), dimension(1:4) :: plane
            type(RectCub_type) :: cell
            call new(cell)
            if (this%vf%cfg%iproc.eq.1) then
               do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
                  do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                     do i=this%vf%cfg%imino,this%vf%cfg%imin-1
                        ! Extract plane data and copy in overlap
                        plane=getPlane(this%vf%liquid_gas_interface(this%vf%cfg%imin,j,k),0)
                        call construct_2pt(cell,[this%vf%cfg%x(i  ),this%vf%cfg%y(j  ),this%vf%cfg%z(k  )],&
                        &                       [this%vf%cfg%x(i+1),this%vf%cfg%y(j+1),this%vf%cfg%z(k+1)])
                        plane(4)=dot_product(plane(1:3),[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)])
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        call matchVolumeFraction(cell,this%vf%VF(i,j,k),this%vf%liquid_gas_interface(i,j,k))
                     end do
                  end do
               end do
            end if
         end block neumann_irl
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         call this%fs%initialize(cfg=this%cfg,name='Two-phase NS')
         this%fs%theta=this%fs%theta+1.0e-2_WP
         ! Set the flow properties
         call param_read('Liquid dynamic viscosity',this%fs%visc_l)
         call param_read('Gas dynamic viscosity'   ,this%fs%visc_g)
         call param_read('Liquid density',this%fs%rho_l)
         call param_read('Gas density'   ,this%fs%rho_g)
         call param_read('Surface tension coefficient',this%fs%sigma)
         ! Inflow on the left
         call this%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp'  ,type=slip           ,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym'  ,type=slip           ,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp'  ,type=slip           ,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm'  ,type=slip           ,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
      end block create_flow_solver
      
      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         use string,     only: str_long
         use messager,   only: log
         character(str_long) :: message
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         real(WP) :: r
         ! Initialize density
         this%resU=this%fs%rho_l*this%vf%VF+this%fs%rho_g*(1.0_WP-this%vf%VF); call this%fs%update_density(rho=this%resU)
         ! Read in inflow conditions
         call param_read('Gas height',this%Hg)
         call param_read('Lip height',this%lip,default=0.0_WP)
         call param_read('Gas velocity',this%Ug)
         call param_read('Liquid velocity',this%Ul)
         ! Compute gas vorticity thickness from Marmottant's correlation
         this%dg=this%Hg*5.6_WP*(this%fs%rho_g*this%Ug*this%Hg/this%fs%visc_g)**(-0.5_WP)
         if (this%fs%cfg%amRoot) then
            write(message,'("[Gas vorticity thickness] => dg =",es12.5)') this%dg; call log(message)
         end if
         ! Apply inflow velocity profile
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            ! Compute radius
            r=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
            ! Set liquid Poiseuille profile
            if (r.le.0.5_WP*this%Dl) this%fs%U(i,j,k)=2.0_WP*this%Ul*(1.0_WP-r/(0.5_WP*this%Dl))**2
            ! Set gas profile
            if (r.ge.0.5_WP*this%Dl+this%lip.and.r.le.0.5_WP*this%Dl+this%lip+this%Hg) this%fs%U(i,j,k)=this%Ug*erf((r-(0.5_WP*this%Dl+this%lip))/this%dg)*erf(((0.5_WP*this%Dl+this%lip+this%Hg)-r)/this%dg)
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Copy to Umid and make it solenoidal
         call this%fs%get_Umid()
         call this%fs%correct_mfr()
         call this%fs%update_laplacian()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%Umid=this%fs%Umid-this%resU/this%fs%sRHOX**2
         this%fs%Vmid=this%fs%Vmid-this%resV/this%fs%sRHOY**2
         this%fs%Wmid=this%fs%Wmid-this%resW/this%fs%sRHOZ**2
         call this%fs%get_U()
         ! Calculate cell-centered velocities and divergence
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',this%use_sgs)
         if (this%use_sgs) this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         character(len=str_medium) :: filename
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         integer :: i,j,k
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call param_read('Restart output period',this%save_evt%tper)
         ! Read in the I/O partition
         call param_read('I/O partition',iopartition)
         ! Check if we are restarting
         call param_read('Restart from',filename,default='')
         this%restarted=.false.; if (len_trim(filename).gt.0) this%restarted=.true.
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata=trim(filename))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                     call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Update the band
            call this%vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Recalculate density
            this%resU=this%fs%rho_l*this%vf%VF+this%fs%rho_g*(1.0_WP-this%vf%VF)
            call this%fs%update_density(rho=this%resU)
            this%fs%sRHOxold=this%fs%sRHOx
            this%fs%sRHOyold=this%fs%sRHOy
            this%fs%sRHOzold=this%fs%sRHOz
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='Umid',var=this%fs%Umid)
            call this%df%pull(name='Vmid',var=this%fs%Vmid)
            call this%df%pull(name='Wmid',var=this%fs%Wmid)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=14)
            this%df%valname=['t ','dt']
            this%df%varname=['U   ','V   ','W   ','Umid','Vmid','Wmid','P   ','Pjx ','Pjy ','Pjz ','P11 ','P12 ','P13 ','P14 ']
         end if
      end block handle_restart
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         this%smesh=surfmesh(nvar=0,name='plic')
         call this%vf%update_surfmesh(this%smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='jet')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_surface('plic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
      
   contains
      
      
      !> Function that defines a level set function for a cylinder
      function levelset_cylinder(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP*this%Dl-sqrt(xyz(2)**2+xyz(3)**2)
      end function levelset_cylinder
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer.or.&
         &   j.le.pg%jmin+this%nlayer.or.&
         &   j.ge.pg%jmax-this%nlayer.or.&
         &   k.le.pg%kmin+this%nlayer.or.&
         &   k.ge.pg%kmax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
      !> Function that localizes the right (x+) of the domain
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
      
      !> Function that localizes the left (x-) of the domain
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator
      
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      
      !> Function that localizes the front (z+) of the domain
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      
      !> Function that localizes the back (z-) of the domain
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      

   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(coaxialjet), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF
      
      ! Remember old velocities and sRHOs
      this%fs%Uold=this%fs%U; this%fs%sRHOxold=this%fs%sRHOx
      this%fs%Vold=this%fs%V; this%fs%sRHOyold=this%fs%sRHOy
      this%fs%Wold=this%fs%W; this%fs%sRHOzold=this%fs%sRHOz
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
            
         ! VOF equation ====================================================
         ! Advance VOF equation
         this%vf%VF=this%vf%VFold
         if (this%time%it.eq.this%time%itmax) then   
            call this%vf%advance(dt=this%time%dt,U=this%fs%Umid,V=this%fs%Vmid,W=this%fs%Wmid)
         else
            call this%vf%advance_tmp(dt=this%time%dt,U=this%fs%Umid,V=this%fs%Vmid,W=this%fs%Wmid)
         end if
         
         ! Update sqrt(face density) and momentum vector
         this%resU=this%fs%rho_l*this%vf%VF+this%fs%rho_g*(1.0_WP-this%vf%VF); call this%fs%update_density(rho=this%resU)
         this%fs%rhoU=this%fs%rho_l*this%vf%UFl(1,:,:,:)+this%fs%rho_g*this%vf%UFg(1,:,:,:)
         this%fs%rhoV=this%fs%rho_l*this%vf%UFl(2,:,:,:)+this%fs%rho_g*this%vf%UFg(2,:,:,:)
         this%fs%rhoW=this%fs%rho_l*this%vf%UFl(3,:,:,:)+this%fs%rho_g*this%vf%UFg(3,:,:,:)
         
         ! Prepare new staggered viscosity (at n+1)
         call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
         
         ! Turbulence modeling
         if (this%use_sgs) then
            sgs_modeling: block
               use sgsmodel_class, only: vreman
               integer :: i,j,k
               this%resU=this%fs%rho_l*this%vf%VF+this%fs%rho_g*(1.0_WP-this%vf%VF)
               call this%fs%get_gradUmid(this%gradU)
               call this%sgs%get_visc(type=vreman,dt=this%time%dt,rho=this%resU,gradu=this%gradU)
               do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_
                  do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_
                     do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
                        this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
                        this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
                        this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
                        this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
                     end do
                  end do
               end do
            end block sgs_modeling
         end if
         
         ! Momentum equation ===============================================
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-(this%fs%U*this%fs%sRHOX**2-this%fs%Uold*this%fs%sRHOXold**2)+this%time%dt*this%resU
         this%resV=-(this%fs%V*this%fs%sRHOY**2-this%fs%Vold*this%fs%sRHOYold**2)+this%time%dt*this%resV
         this%resW=-(this%fs%W*this%fs%sRHOZ**2-this%fs%Wold*this%fs%sRHOZold**2)+this%time%dt*this%resW
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Compute predictor U
         this%fs%U=this%fs%U+this%resU
         this%fs%V=this%fs%V+this%resV
         this%fs%W=this%fs%W+this%resW
         
         ! Sync and apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Poisson equation ================================================
         ! Compute Umid from U and Uold
         call this%fs%get_Umid()
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct pressure and Umid
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%Umid=this%fs%Umid-this%time%dt*this%resU/((this%fs%sRHOX+this%fs%sRHOXold*(1.0_WP-this%fs%theta)/this%fs%theta)*this%fs%sRHOX)
         this%fs%Vmid=this%fs%Vmid-this%time%dt*this%resV/((this%fs%sRHOY+this%fs%sRHOYold*(1.0_WP-this%fs%theta)/this%fs%theta)*this%fs%sRHOY)
         this%fs%Wmid=this%fs%Wmid-this%time%dt*this%resW/((this%fs%sRHOZ+this%fs%sRHOZold*(1.0_WP-this%fs%theta)/this%fs%theta)*this%fs%sRHOZ)
         
         ! Regenerate U from Umid and Uold
         call this%fs%get_U()
         
         ! Increment sub-iteration counter =================================
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         ! Adjust VOF field
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            if (n.le.this%vof_removal_layer%n_) this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
         ! Also adjust density
         this%resU=this%fs%rho_l*this%vf%VF+this%fs%rho_g*(1.0_WP-this%vf%VF); call this%fs%update_density(rho=this%resU)
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         call this%vf%update_surfmesh(this%smesh)
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'   ,val=this%time%t )
            call this%df%push(name='dt'  ,val=this%time%dt)
            call this%df%push(name='U'   ,var=this%fs%U   )
            call this%df%push(name='V'   ,var=this%fs%V   )
            call this%df%push(name='W'   ,var=this%fs%W   )
            call this%df%push(name='Umid',var=this%fs%Umid)
            call this%df%push(name='Vmid',var=this%fs%Vmid)
            call this%df%push(name='Wmid',var=this%fs%Wmid)
            call this%df%push(name='P'   ,var=this%fs%P   )
            call this%df%push(name='Pjx' ,var=this%fs%Pjx )
            call this%df%push(name='Pjy' ,var=this%fs%Pjy )
            call this%df%push(name='Pjz' ,var=this%fs%Pjz )
            call this%df%push(name='P11' ,var=P11         )
            call this%df%push(name='P12' ,var=P12         )
            call this%df%push(name='P13' ,var=P13         )
            call this%df%push(name='P14' ,var=P14         )
            call this%df%write(fdata='restart/jet_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14)
         end block save_restart
      end if
      
   end subroutine step
   
   
   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(coaxialjet), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      
   end subroutine final
   
   
end module coaxialjet_class