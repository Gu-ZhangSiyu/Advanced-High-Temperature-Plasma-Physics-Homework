!***************************************************************
!	Coded By Siyu Zhang     May. 22, 2026
!
!   Full-PIC 主程序 高温等离子体物理学课题用
!
!	基于 Full-PIC_v2.2.3 版本修改 (Siyu Zhang, April. 22, 2026)
!*********************************************************************
module mesh_parameter
    implicit none
	integer, parameter :: nx = 257, ny = 257, nz = 3		! 网格数
    real(8), parameter :: x_min = -6.4d-2		! 物理边界
    real(8), parameter :: x_max =  6.4d-2
    real(8), parameter :: y_min = -6.4d-2
    real(8), parameter :: y_max =  6.4d-2
    real(8), parameter :: z_min = -6.4d-2
    real(8), parameter :: z_max =  6.4d-2
    real(8), parameter :: dx = (x_max - x_min) / dble(nx - 1)	! 网格间距(均匀网格)
    real(8), parameter :: dy = (y_max - y_min) / dble(ny - 1)
    real(8), parameter :: dz = (z_max - z_min) / dble(nz - 1)	
end module mesh_parameter
!*********************************************************************
module physicals_parameter
    implicit none
	real(8), parameter :: pi = 3.141592653589793d0	! 圆周率
	real(8), parameter :: q  = 1.60217733d-19       ! 元电荷 (C)
	real(8), parameter :: c  = 2.99792458d8			! 光速
	real(8), parameter :: eps0 = 8.8541878128d-12   ! 真空介电常数
	real(8), parameter :: M_ELEC = 9.1093837139d-31	! 电子质量
	real(8), parameter :: M_ION = 25.0d0 * M_ELEC	! 离子质量
end module physicals_parameter
!*********************************************************************
module pic_types
    use physicals_parameter
    implicit none
	type :: particle
		integer(kind=8) :: id	! 粒子ID
		real(8) :: x, y, z		! 粒子位置 (m)
		real(8) :: vx, vy, vz	! 粒子速度 (m/s)
		real(8) :: weight		! 超粒子权重
	end type particle
end module pic_types
!*********************************************************************
module control_parameter
    implicit none
    integer(kind=8) :: nIon  = 50000000	! 初始超离子数
    integer(kind=8) :: nElec = 50000000	! 初始超电子数
    integer :: nStepsDyn = 20000			! FDTD蛙跳法主循环步数
    real(8) :: dt = 0.5d-12					! 时间步长(s)
    integer :: freq_field = 100				! 输出场频率
end module control_parameter
!***************************************************************
module diagnostics_solver
    use mpi; use mesh_parameter; use physicals_parameter; use pic_types
    implicit none
contains
	subroutine Output_Controller(idyn, freq_field, Bx, By)
        implicit none
        integer, intent(in) :: idyn, freq_field
        real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:)
        if (mod(idyn, freq_field) == 0 .or. idyn == 1) then
            call dump_electromagnetic_fields(Bx, By, idyn)
        end if
    end subroutine Output_Controller
!*********************************************************************
subroutine dump_electromagnetic_fields(Bx, By, step)    ! 场输出
    use mesh_parameter, only: nx, ny, nz; use mpi
    implicit none
    integer, intent(in) :: step
    integer :: rank, ierr, u_field, mid_z, i, j
    character(len=64) :: filename
	real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:)
	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
	mid_z = (nz + 1) / 2
	if (rank == 0) then
        write(filename, '("dynamic_fields_Bx/", I0, ".txt")') step
        open(newunit=u_field, file=filename, status='replace')
        do j = 1, ny; write(u_field, '(*(E15.6))') (Bx(i, j, mid_z), i = 1, nx); end do; close(u_field)
        write(filename, '("dynamic_fields_By/", I0, ".txt")') step
        open(newunit=u_field, file=filename, status='replace')
        do j = 1, ny; write(u_field, '(*(E15.6))') (By(i, j, mid_z), i = 1, nx); end do; close(u_field)
    end if
end subroutine dump_electromagnetic_fields
end module diagnostics_solver
!***************************************************************
module data_loader_module
	use mesh_parameter, only: nx, ny, nz
    use pic_types, only: particle
    use mpi
    implicit none
    private
    public :: load_equilibrium, read_particle_data
contains
!*********************************************************************
subroutine load_equilibrium(Bx, By, Bz, Ex, Ey, Ez)		! 读背景场
    implicit none
    real(8), intent(out) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:), Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
    integer :: io_stat
    real(8), allocatable :: temp_3d(:,:,:)
    allocate(temp_3d(nx, ny, nz))
    open(unit=10, file='out_3d/Bx_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) temp_3d; close(10); Bx(1:nx, 1:ny, 1:nz) = temp_3d
    open(unit=10, file='out_3d/By_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) temp_3d; close(10); By(1:nx, 1:ny, 1:nz) = temp_3d
    open(unit=10, file='out_3d/Bz_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) temp_3d; close(10); Bz(1:nx, 1:ny, 1:nz) = temp_3d
    open(unit=10, file='out_3d/Ex_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) temp_3d; close(10); Ex(1:nx, 1:ny, 1:nz) = temp_3d
    open(unit=10, file='out_3d/Ey_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) temp_3d; close(10); Ey(1:nx, 1:ny, 1:nz) = temp_3d
    open(unit=10, file='out_3d/Ez_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) temp_3d; close(10); Ez(1:nx, 1:ny, 1:nz) = temp_3d
    deallocate(temp_3d)
end subroutine load_equilibrium
!*********************************************************************
subroutine read_particle_data(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, ierr)	! 读粒子
    implicit none
    type(particle), intent(out) :: ion_local(:), elec_local(:)
    integer, intent(in) :: nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, ierr
    character(len=64) :: filename
    write(filename, '(A,A)') 'init_particles/ion_x', '.bin'
	call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'x', rank, nprocs, ierr)
    write(filename, '(A,A)') 'init_particles/ion_y', '.bin'
    call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'y', rank, nprocs, ierr)
	write(filename, '(A,A)') 'init_particles/ion_z', '.bin'
    call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'z', rank, nprocs, ierr)
    write(filename, '(A,A)') 'init_particles/ion_vx', '.bin'
    call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'vx', rank, nprocs, ierr)
    write(filename, '(A,A)') 'init_particles/ion_vy', '.bin'
    call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'vy', rank, nprocs, ierr)   
    write(filename, '(A,A)') 'init_particles/ion_vz', '.bin'
    call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'vz', rank, nprocs, ierr)   
    write(filename, '(A,A)') 'init_particles/ion_weight', '.bin'
    call read_component_parallel(ion_local, nIon_local, ion_offset, filename, 'weight', rank, nprocs, ierr)
    write(filename, '(A,A)') 'init_particles/elec_x', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'x', rank, nprocs, ierr)  
	write(filename, '(A,A)') 'init_particles/elec_y', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'y', rank, nprocs, ierr) 
    write(filename, '(A,A)') 'init_particles/elec_z', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'z', rank, nprocs, ierr)    
    write(filename, '(A,A)') 'init_particles/elec_vx', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'vx', rank, nprocs, ierr)   
    write(filename, '(A,A)') 'init_particles/elec_vy', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'vy', rank, nprocs, ierr)    
    write(filename, '(A,A)') 'init_particles/elec_vz', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'vz', rank, nprocs, ierr)    
    write(filename, '(A,A)') 'init_particles/elec_weight', '.bin'
    call read_component_parallel(elec_local, nElec_local, elec_offset, filename, 'weight', rank, nprocs, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, 'Particle loading complete.'
end subroutine read_particle_data
!*********************************************************************
subroutine read_component_parallel(part, nLocal, offset, filename, component, rank, nprocs, ierr)	! 读粒子辅助
    implicit none
    type(particle), intent(inout) :: part(:)	! 进程本地粒子数组
    integer, intent(in) :: nLocal, offset, rank, nprocs, ierr
    character(len=*), intent(in) :: filename, component
    real(8), allocatable :: local_data(:)
    integer :: p, i
    integer, allocatable :: all_nLocal(:)
    integer(kind=8), allocatable :: byte_offsets(:)
    
    allocate(local_data(nLocal))
    allocate(all_nLocal(0:nprocs-1))
    allocate(byte_offsets(0:nprocs-1))
    
    call MPI_Allgather(nLocal, 1, MPI_INTEGER, all_nLocal, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    
	byte_offsets(0) = 0_8	! 强制指定双精度浮点数为 8 字节
    do p = 1, nprocs-1
        byte_offsets(p) = byte_offsets(p-1) + int(all_nLocal(p-1), 8) * 8_8
    end do
    
    ! 每个进程打开文件寻址数据块
	open(unit=30, file=filename, status='old', access='stream', action='read', form='unformatted')
        read(30, pos=byte_offsets(rank) + 1_8) local_data
    close(30)
    deallocate(all_nLocal, byte_offsets)
    
	select case (trim(component))
		case ('x'); do i = 1, nLocal; part(i)%x = local_data(i); end do
		case ('y'); do i = 1, nLocal; part(i)%y = local_data(i); end do
		case ('z'); do i = 1, nLocal; part(i)%z = local_data(i); end do
		case ('vx'); do i = 1, nLocal; part(i)%vx = local_data(i); end do
		case ('vy'); do i = 1, nLocal; part(i)%vy = local_data(i); end do
		case ('vz'); do i = 1, nLocal; part(i)%vz = local_data(i); end do
		case ('weight'); do i = 1, nLocal; part(i)%weight = local_data(i); end do
		case default	
			if (rank == 0) print *, 'ERROR: Unknown component: ', trim(component)
			call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end select
    deallocate(local_data)
end subroutine read_component_parallel
end module data_loader_module
!*********************************************************************
program pic_equilibrium
    use mpi; use iso_c_binding, only: c_ptr, c_f_pointer
	use mesh_parameter; use physicals_parameter; use pic_types; use control_parameter
    use particle_solver_yee, only: Push_Particle
    use diagnostics_solver, only: Output_Controller
	use data_loader_module, only: load_equilibrium, read_particle_data
    implicit none
	integer :: rank, nprocs, ierr, node_comm, node_rank, disp_unit, win_Bx, win_By, win_Bz, win_Ex, win_Ey, win_Ez, &
				nIon_local, nElec_local, ion_offset, elec_offset, base_Ion, extra_Ion, base_Elec, extra_Elec, &
				i, j, k, idyn
	integer(kind=8) :: global_nIon_alive, global_nElec_alive
	integer(kind=MPI_ADDRESS_KIND) :: window_size
	real(8), pointer :: Bx0(:,:,:), By0(:,:,:), Bz0(:,:,:), Ex0(:,:,:), Ey0(:,:,:), Ez0(:,:,:), temp_ptr(:,:,:)
    real(8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:), Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), Jx(:,:,:), Jy(:,:,:), Jz(:,:,:), &
							Bx_half(:,:,:), By_half(:,:,:), Bz_half(:,:,:), Bx_old(:,:,:), By_old(:,:,:), Bz_old(:,:,:)
	real(8) :: qm_ion, qm_elec
    type(c_ptr) :: base_Bx, base_By, base_Bz, base_Ex, base_Ey, base_Ez
    type(particle), pointer :: ion_local(:), elec_local(:)

	! 初始化MPI
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierr)
    call MPI_Comm_rank(node_comm, node_rank, ierr)
	if (rank == 0) print *, '----------------- Start Code -----------------'
	
	! 申请物理内存
    if (node_rank == 0) then
        window_size = int(nx+2, MPI_ADDRESS_KIND) * int(ny+2, MPI_ADDRESS_KIND) * int(nz+2, MPI_ADDRESS_KIND) * 8_MPI_ADDRESS_KIND
    else
        window_size = 0_MPI_ADDRESS_KIND
    end if
    disp_unit = 8

    ! 分配背景场并进行指针重映射
    call MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, node_comm, base_Bx, win_Bx, ierr)
    if (node_rank /= 0) call MPI_Win_shared_query(win_Bx, 0, window_size, disp_unit, base_Bx, ierr)
    call c_f_pointer(base_Bx, temp_ptr, [nx+2, ny+2, nz+2])
    Bx0(0:nx+1, 0:ny+1, 0:nz+1) => temp_ptr
    call MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, node_comm, base_By, win_By, ierr)
    if (node_rank /= 0) call MPI_Win_shared_query(win_By, 0, window_size, disp_unit, base_By, ierr)
    call c_f_pointer(base_By, temp_ptr, [nx+2, ny+2, nz+2])
    By0(0:nx+1, 0:ny+1, 0:nz+1) => temp_ptr
    call MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, node_comm, base_Bz, win_Bz, ierr)
    if (node_rank /= 0) call MPI_Win_shared_query(win_Bz, 0, window_size, disp_unit, base_Bz, ierr)
    call c_f_pointer(base_Bz, temp_ptr, [nx+2, ny+2, nz+2])
    Bz0(0:nx+1, 0:ny+1, 0:nz+1) => temp_ptr
    call MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, node_comm, base_Ex, win_Ex, ierr)
    if (node_rank /= 0) call MPI_Win_shared_query(win_Ex, 0, window_size, disp_unit, base_Ex, ierr)
    call c_f_pointer(base_Ex, temp_ptr, [nx+2, ny+2, nz+2])
    Ex0(0:nx+1, 0:ny+1, 0:nz+1) => temp_ptr
    call MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, node_comm, base_Ey, win_Ey, ierr)
    if (node_rank /= 0) call MPI_Win_shared_query(win_Ey, 0, window_size, disp_unit, base_Ey, ierr)
    call c_f_pointer(base_Ey, temp_ptr, [nx+2, ny+2, nz+2])
    Ey0(0:nx+1, 0:ny+1, 0:nz+1) => temp_ptr
    call MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, node_comm, base_Ez, win_Ez, ierr)
    if (node_rank /= 0) call MPI_Win_shared_query(win_Ez, 0, window_size, disp_unit, base_Ez, ierr)
    call c_f_pointer(base_Ez, temp_ptr, [nx+2, ny+2, nz+2])
    Ez0(0:nx+1, 0:ny+1, 0:nz+1) => temp_ptr

    if (node_rank == 0) then
        Bx0 = 0.0d0; By0 = 0.0d0; Bz0 = 0.0d0; Ex0 = 0.0d0; Ey0 = 0.0d0; Ez0 = 0.0d0
    end if
    call MPI_Barrier(node_comm, ierr)
	
	! 分配动态场数组(包含卫星层)
    allocate(Bx(0:nx+1, 0:ny+1, 0:nz+1), By(0:nx+1, 0:ny+1, 0:nz+1), Bz(0:nx+1, 0:ny+1, 0:nz+1), &
            Ex(0:nx+1, 0:ny+1, 0:nz+1), Ey(0:nx+1, 0:ny+1, 0:nz+1), Ez(0:nx+1, 0:ny+1, 0:nz+1), &
            Bx_half(0:nx+1, 0:ny+1, 0:nz+1), By_half(0:nx+1, 0:ny+1, 0:nz+1), Bz_half(0:nx+1, 0:ny+1, 0:nz+1), &
            Bx_old(0:nx+1, 0:ny+1, 0:nz+1), By_old(0:nx+1, 0:ny+1, 0:nz+1), Bz_old(0:nx+1, 0:ny+1, 0:nz+1), &
            Jx(0:nx+1, 0:ny+1, 0:nz+1), Jy(0:nx+1, 0:ny+1, 0:nz+1), Jz(0:nx+1, 0:ny+1, 0:nz+1))
	
    Bx = 0.0d0; By = 0.0d0; Bz = 0.0d0; Ex = 0.0d0; Ey = 0.0d0; Ez = 0.0d0; Jx = 0.0d0; Jy = 0.0d0; Jz = 0.0d0
    Bx_half = 0.0d0; By_half = 0.0d0; Bz_half = 0.0d0; Bx_old = 0.0d0; By_old = 0.0d0; Bz_old = 0.0d0

    ! 重构粒子数
	global_nIon_alive = nIon; global_nElec_alive = nElec

	! 计算本地粒子数和偏移量
	base_Ion = global_nIon_alive / nprocs
    extra_Ion = mod(global_nIon_alive, nprocs)
	if (rank < extra_Ion) then
		nIon_local = base_Ion + 1
		ion_offset = rank * (base_Ion + 1)
	else
		nIon_local = base_Ion
		ion_offset = extra_Ion * (base_Ion + 1) + (rank - extra_Ion) * base_Ion
	end if
	base_Elec = global_nElec_alive / nprocs
    extra_Elec = mod(global_nElec_alive, nprocs)
	if (rank < extra_Elec) then
		nElec_local = base_Elec + 1
		elec_offset = rank * (base_Elec + 1)
	else
		nElec_local = base_Elec
		elec_offset = extra_Elec * (base_Elec + 1) + (rank - extra_Elec) * base_Elec
	end if
	
	! 分配本地粒子数组
	allocate(ion_local(nIon_local))
	allocate(elec_local(nElec_local))

	! 加载初始场并广播
    if (node_rank == 0) then
        call load_equilibrium(Bx0, By0, Bz0, Ex0, Ey0, Ez0)
        call Ghost_Cells(Bx0, 'Bx'); call Ghost_Cells(By0, 'By'); call Ghost_Cells(Bz0, 'Bz')
        call Ghost_Cells(Ex0, 'Ex'); call Ghost_Cells(Ey0, 'Ey'); call Ghost_Cells(Ez0, 'Ez')
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
	
    ! 初始化总场
	Bx = Bx0; By = By0; Bz = Bz0; Ex = Ex0; Ey = Ey0; Ez = Ez0;

    ! 初始化卫星层
    call Ghost_Cells(Bx, 'Bx'); call Ghost_Cells(By, 'By'); call Ghost_Cells(Bz, 'Bz')
    call Ghost_Cells(Ex, 'Ex'); call Ghost_Cells(Ey, 'Ey'); call Ghost_Cells(Ez, 'Ez')

	! 加载粒子
	call read_particle_data(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, ierr)
	! 定义粒子ID
	do i = 1, nIon_local
		ion_local(i)%id = int(ion_offset, 8) + int(i, 8)
	end do
	do i = 1, nElec_local
		elec_local(i)%id = int(elec_offset, 8) + int(i, 8)
	end do

    qm_ion = q / M_ION; qm_elec = -q / M_ELEC   ! 预计算荷质比

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (rank == 0) print *, '----------------- Start Main -----------------'

    ! 磁场时间与读取初始值同步
    call Update_B_FDTD(Bx, By, Bz, Ex, Ey, Ez, -dt/2.0d0)
    call Ghost_Cells(Bx, 'Bx'); call Ghost_Cells(By, 'By'); call Ghost_Cells(Bz, 'Bz')

    ! 蛙跳主循环
    do idyn = 1, nStepsDyn
        if (rank == 0 .and. mod(idyn, 100) == 0) print *, 'Step: ', idyn

        ! 备份旧磁场 B(n-1/2)
        Bx_old = Bx; By_old = By; Bz_old = Bz

        ! 更新磁场 B(n-1/2) -> B(n+1/2)
        call Update_B_FDTD(Bx, By, Bz, Ex, Ey, Ez, dt)
        call Ghost_Cells(Bx, 'Bx'); call Ghost_Cells(By, 'By'); call Ghost_Cells(Bz, 'Bz')

        ! 构造时间中心磁场 B(n)
        Bx_half = 0.5d0 * (Bx + Bx_old); By_half = 0.5d0 * (By + By_old); Bz_half = 0.5d0 * (Bz + Bz_old)

        ! 粒子推进 x(n), v(n-1/2) -> x(n+1), v(n+1/2)
		call Run_Particle_Pusher(ion_local, nIon_local, Ex, Ey, Ez, Bx_half, By_half, Bz_half, dt, 0, qm_ion)
        call Run_Particle_Pusher(elec_local, nElec_local, Ex, Ey, Ez, Bx_half, By_half, Bz_half, dt, 1, qm_elec)

        ! 电流沉积 J(n+1/2)
        call Deposit_Currents(ion_local, nIon_local, elec_local, nElec_local, Jx, Jy, Jz, dt)

        ! 空间平滑滤波, 扣除本底电流
        call Filter_Currents_3D(Jx, Jy, Jz)

        ! 更新电场 E(n) -> E(n+1)
        call Update_E_FDTD(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, dt)
        call Ghost_Cells(Ex, 'Ex'); call Ghost_Cells(Ey, 'Ey'); call Ghost_Cells(Ez, 'Ez')

        ! 输出				
		call Output_Controller(idyn=idyn, freq_field=freq_field, Bx=Bx, By=By)
    end do
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (rank == 0) print *, '----------------- End Main -----------------'
	call MPI_Finalize(ierr)
contains
!*********************************************************************
subroutine Run_Particle_Pusher(part_array, nP, Ex, Ey, Ez, Bx, By, Bz, dt, rel_switch, qm_ratio_species)	! 调用粒子推进求解器
    use mesh_parameter; use pic_types
    implicit none
    type(particle), intent(inout) :: part_array(:)
    real(8), intent(in) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:), Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:)
    real(8), intent(in) :: dt, qm_ratio_species
    integer, intent(in) :: nP, rel_switch	! 0=Non-Rel, 1=Rel
    integer :: i, status_code
    real(8) :: p_in(3), v_in(3), p_out(3), v_out(3), origin_vec(3), step_vec(3)
    real(8) :: qm_ratio
    integer :: dim_vec(3)

    origin_vec = [x_min, y_min, z_min]
    step_vec = [dx, dy, dz]
    dim_vec = [nx, ny, nz]

    do i = 1, nP
        if (part_array(i)%weight <= 0.0d0) cycle

        ! 提取数据到临时数组
        p_in = [part_array(i)%x, part_array(i)%y, part_array(i)%z]
        v_in = [part_array(i)%vx, part_array(i)%vy, part_array(i)%vz]

        ! 调用boris求解器
        call Push_Particle(pos_in=p_in, vel_in=v_in, qm=qm_ratio_species, dt=dt, Ex=Ex, Ey=Ey, Ez=Ez, Bx=Bx, By=By, Bz=Bz, &
                            grid_origin=origin_vec, grid_steps=step_vec, grid_dims=dim_vec, pos_out=p_out, vel_out=v_out, &
                            status=status_code, boundary_control=0, relativistic_control=rel_switch)

        if (status_code /= 0) then
            part_array(i)%weight = 0.0d0
        else
			if (p_out(1) >= x_max) p_out(1) = p_out(1) - (x_max - x_min)
            if (p_out(1) < x_min)  p_out(1) = p_out(1) + (x_max - x_min)
            if (p_out(3) >= z_max) p_out(3) = p_out(3) - (z_max - z_min)
            if (p_out(3) < z_min)  p_out(3) = p_out(3) + (z_max - z_min)
            if (p_out(2) >= y_max .or. p_out(2) <= y_min) part_array(i)%weight = 0.0d0
            part_array(i)%x  = p_out(1)
            part_array(i)%y  = p_out(2)
            part_array(i)%z  = p_out(3)
            part_array(i)%vx = v_out(1)
            part_array(i)%vy = v_out(2)
            part_array(i)%vz = v_out(3)
        end if
    end do
end subroutine Run_Particle_Pusher
!*********************************************************************
subroutine Deposit_Currents(ions, nIon_local, electrons, nElec_local, Jx, Jy, Jz, dt)   ! 电流沉积
    use mpi; use mesh_parameter; use physicals_parameter; use pic_types
    implicit none
    type(particle), intent(in) :: ions(:), electrons(:)
    integer, intent(in) :: nIon_local, nElec_local
    real(8), intent(in) :: dt
    real(8), intent(out) :: Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    real(8), allocatable, save :: J_work(:, :, :)
    real(8) :: node_vol_inv
    integer :: ierr
    
    Jx = 0.0d0; Jy = 0.0d0; Jz = 0.0d0
    node_vol_inv = 1.0d0 / (dx * dy * dz)

    call Deposit_Species_Yee(ions, nIon_local, Jx, Jy, Jz, q, dt)
    call Deposit_Species_Yee(electrons, nElec_local, Jx, Jy, Jz, -q, dt)

    if (.not. allocated(J_work)) allocate(J_work(nx, ny, nz))
    J_work = Jx(1:nx, 1:ny, 1:nz) ! 提取
    call MPI_Allreduce(MPI_IN_PLACE, J_work, nx*ny*nz, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Jx(1:nx, 1:ny, 1:nz) = J_work ! 填回
    J_work = Jy(1:nx, 1:ny, 1:nz)
    call MPI_Allreduce(MPI_IN_PLACE, J_work, nx*ny*nz, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Jy(1:nx, 1:ny, 1:nz) = J_work   
    J_work = Jz(1:nx, 1:ny, 1:nz)
    call MPI_Allreduce(MPI_IN_PLACE, J_work, nx*ny*nz, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Jz(1:nx, 1:ny, 1:nz) = J_work   
    Jx = Jx * node_vol_inv; Jy = Jy * node_vol_inv; Jz = Jz * node_vol_inv  ! 归一化

    Jx(1,:,:) = Jx(1,:,:) + Jx(nx,:,:)
    Jx(nx-1,:,:) = Jx(nx-1,:,:) + Jx(0,:,:)
    Jx(2,:,:) = Jx(2,:,:) + Jx(nx+1,:,:)    
    Jy(1,:,:) = Jy(1,:,:) + Jy(nx,:,:)
    Jy(nx-1,:,:) = Jy(nx-1,:,:) + Jy(0,:,:)
    Jy(2,:,:) = Jy(2,:,:) + Jy(nx+1,:,:)    
    Jz(1,:,:) = Jz(1,:,:) + Jz(nx,:,:)
    Jz(nx-1,:,:) = Jz(nx-1,:,:) + Jz(0,:,:)
    Jz(2,:,:) = Jz(2,:,:) + Jz(nx+1,:,:)

    if (nz == 3) then
        Jx(:,:,1) = Jx(:,:,1) + Jx(:,:,3)
        Jy(:,:,1) = Jy(:,:,1) + Jy(:,:,3)
        Jz(:,:,1) = Jz(:,:,1) + Jz(:,:,3)        
        Jx(:,:,2) = Jx(:,:,2) + Jx(:,:,0) + Jx(:,:,4)
        Jy(:,:,2) = Jy(:,:,2) + Jy(:,:,0) + Jy(:,:,4)
        Jz(:,:,2) = Jz(:,:,2) + Jz(:,:,0) + Jz(:,:,4)
    else
        Jx(:,:,1) = Jx(:,:,1) + Jx(:,:,nz)
        Jx(:,:,nz-1) = Jx(:,:,nz-1) + Jx(:,:,0)
        Jx(:,:,2) = Jx(:,:,2) + Jx(:,:,nz+1)
        Jy(:,:,1) = Jy(:,:,1) + Jy(:,:,nz)
        Jy(:,:,nz-1) = Jy(:,:,nz-1) + Jy(:,:,0)
        Jy(:,:,2) = Jy(:,:,2) + Jy(:,:,nz+1)        
        Jz(:,:,1) = Jz(:,:,1) + Jz(:,:,nz)
        Jz(:,:,nz-1) = Jz(:,:,nz-1) + Jz(:,:,0)
        Jz(:,:,2) = Jz(:,:,2) + Jz(:,:,nz+1)
    end if
    call Ghost_Cells(Jx, 'Jx'); call Ghost_Cells(Jy, 'Jy'); call Ghost_Cells(Jz, 'Jz')
end subroutine Deposit_Currents
!*********************************************************************
subroutine Deposit_Species_Yee(part, nP, Jx_out, Jy_out, Jz_out, species_charge, dt)    ! 保电荷VB沉积
    use mesh_parameter; use pic_types
    implicit none
    type(particle), intent(in) :: part(:)
    integer, intent(in) :: nP
    real(8), intent(inout) :: Jx_out(0:,0:,0:), Jy_out(0:,0:,0:), Jz_out(0:,0:,0:)
    real(8), intent(in) :: species_charge, dt
    integer :: ip, iter, iu, iv, iw
    real(8) :: q_w, dx_inv, dy_inv, dz_inv, tu, tv, tw, t_min, u1, v1, w1, u2, v2, w2
    real(8) :: u_curr, v_curr, w_curr, u_next, v_next, w_next, u_min, u_max, v_min, v_max, w_min, w_max
    real(8), parameter :: eps = 1.0d-12

    dx_inv = 1.0d0 / dx; dy_inv = 1.0d0 / dy; dz_inv = 1.0d0 / dz
    do ip = 1, nP
        if (part(ip)%weight <= 0.0d0) cycle
        q_w = species_charge * part(ip)%weight

        ! 转换到网格归一化坐标系
        u2 = (part(ip)%x - x_min) * dx_inv
        v2 = (part(ip)%y - y_min) * dy_inv
        w2 = (part(ip)%z - z_min) * dz_inv

        ! v(n+1/2)回溯出n时刻的位置
        u1 = u2 - part(ip)%vx * dt * dx_inv
        v1 = v2 - part(ip)%vy * dt * dy_inv
        w1 = w2 - part(ip)%vz * dt * dz_inv

        if (u2 < 0.0d0 .or. u2 >= dble(nx) .or. v2 < 0.0d0 .or. v2 >= dble(ny) .or. w2 < 0.0d0 .or. w2 >= dble(nz)) cycle
        u_curr = u1; v_curr = v1; w_curr = w1

        ! 射线步进算法
        do iter = 1, 4
            ! 确定当前线段所在的基准网格系
            iu = floor(u_curr + eps)
            iv = floor(v_curr + eps)
            iw = floor(w_curr + eps)

            ! 修正落在边界且速度向负方向的边界判定
            if (abs(u_curr - dble(iu)) < eps .and. u2 < u_curr) iu = iu - 1
            if (abs(v_curr - dble(iv)) < eps .and. v2 < v_curr) iv = iv - 1
            if (abs(w_curr - dble(iw)) < eps .and. w2 < w_curr) iw = iw - 1
            u_min = dble(iu); u_max = dble(iu + 1)
            v_min = dble(iv); v_max = dble(iv + 1)
            w_min = dble(iw); w_max = dble(iw + 1)

            ! 计算与各个胞面相交所需的参数步长t 
            tu = 1.0d0; tv = 1.0d0; tw = 1.0d0
            if (u2 > u_max) tu = (u_max - u_curr) / (u2 - u_curr)
            if (u2 < u_min) tu = (u_min - u_curr) / (u2 - u_curr)
            if (v2 > v_max) tv = (v_max - v_curr) / (v2 - v_curr)
            if (v2 < v_min) tv = (v_min - v_curr) / (v2 - v_curr)
            if (w2 > w_max) tw = (w_max - w_curr) / (w2 - w_curr)
            if (w2 < w_min) tw = (w_min - w_curr) / (w2 - w_curr)
            t_min = min(tu, tv, tw, 1.0d0)

            ! 当前胞元内的局部截断段
            u_next = u_curr + t_min * (u2 - u_curr)
            v_next = v_curr + t_min * (v2 - v_curr)
            w_next = w_curr + t_min * (w2 - w_curr)

            call Deposit_Segment_VB(u_curr, v_curr, w_curr, u_next, v_next, w_next, iu, iv, iw, q_w, dt, Jx_out, Jy_out, Jz_out)
            if (t_min >= 1.0d0) exit
            u_curr = u_next; v_curr = v_next; w_curr = w_next
        end do
    end do
end subroutine Deposit_Species_Yee
!*********************************************************************
subroutine Deposit_Segment_VB(uA, vA, wA, uB, vB, wB, iu, iv, iw, q_w, dt, Jx_out, Jy_out, Jz_out)  !单胞几何体积通量计算(Villasenor-Buneman1992)
    use mesh_parameter, only: dx, dy, dz, nx, ny, nz
    implicit none
    real(8), intent(in) :: uA, vA, wA, uB, vB, wB, q_w, dt
    integer, intent(in) :: iu, iv, iw
    real(8), intent(inout) :: Jx_out(0:,0:,0:), Jy_out(0:,0:,0:), Jz_out(0:,0:,0:)
    real(8) :: du, dv, dw, u_mid, v_mid, w_mid, Wx0, Wx1, Wy0, Wy1, Wz0, Wz1
    real(8) :: Fx00, Fx01, Fx10, Fx11, Fy00, Fy01, Fy10, Fy11, Fz00, Fz01, Fz10, Fz11
    integer :: i0, j0, k0

    du = uB - uA; dv = vB - vA; dw = wB - wA    ! 线段空间位移量(归一化)

    ! 线段中心点在该胞内的相对局部坐标 [0, 1)
    u_mid = 0.5d0 * (uA + uB) - dble(iu)
    v_mid = 0.5d0 * (vA + vB) - dble(iv)
    w_mid = 0.5d0 * (wA + wB) - dble(iw)

    ! 映射至Fortran数组索引
    i0 = iu + 1; j0 = iv + 1; k0 = iw + 1

    if (i0 < 0 .or. i0 > nx .or. j0 < 1 .or. j0 >= ny .or. k0 < 0 .or. k0 > nz) return

    ! 横向分布形函数权重 (一阶CIC多项式解析解)
    Wx1 = u_mid; Wx0 = 1.0d0 - Wx1
    Wy1 = v_mid; Wy0 = 1.0d0 - Wy1
    Wz1 = w_mid; Wz0 = 1.0d0 - Wz1

    ! 构造精确体积通量(Fx, Fy, Fz)
    Fx00 = du * (Wy0 * Wz0 + (dv * dw) / 12.0d0)
    Fx01 = du * (Wy0 * Wz1 - (dv * dw) / 12.0d0)
    Fx10 = du * (Wy1 * Wz0 - (dv * dw) / 12.0d0)
    Fx11 = du * (Wy1 * Wz1 + (dv * dw) / 12.0d0)
    Fy00 = dv * (Wx0 * Wz0 + (du * dw) / 12.0d0)
    Fy01 = dv * (Wx0 * Wz1 - (du * dw) / 12.0d0)
    Fy10 = dv * (Wx1 * Wz0 - (du * dw) / 12.0d0)
    Fy11 = dv * (Wx1 * Wz1 + (du * dw) / 12.0d0)
    Fz00 = dw * (Wx0 * Wy0 + (du * dv) / 12.0d0)
    Fz01 = dw * (Wx0 * Wy1 - (du * dv) / 12.0d0)
    Fz10 = dw * (Wx1 * Wy0 - (du * dv) / 12.0d0)
    Fz11 = dw * (Wx1 * Wy1 + (du * dv) / 12.0d0)

    ! 叠加到全局电流场
    Jx_out(i0, j0, k0) = Jx_out(i0, j0, k0) + q_w * (dx / dt) * Fx00
    Jx_out(i0, j0, k0+1) = Jx_out(i0, j0, k0+1) + q_w * (dx / dt) * Fx01
    Jx_out(i0, j0+1, k0) = Jx_out(i0, j0+1, k0) + q_w * (dx / dt) * Fx10
    Jx_out(i0, j0+1, k0+1) = Jx_out(i0, j0+1, k0+1) + q_w * (dx / dt) * Fx11
    Jy_out(i0, j0, k0) = Jy_out(i0, j0, k0) + q_w * (dy / dt) * Fy00
    Jy_out(i0, j0, k0+1) = Jy_out(i0, j0, k0+1) + q_w * (dy / dt) * Fy01
    Jy_out(i0+1, j0, k0) = Jy_out(i0+1, j0, k0) + q_w * (dy / dt) * Fy10
    Jy_out(i0+1, j0, k0+1) = Jy_out(i0+1, j0, k0+1) + q_w * (dy / dt) * Fy11
    Jz_out(i0, j0, k0) = Jz_out(i0, j0, k0) + q_w * (dz / dt) * Fz00
    Jz_out(i0, j0+1, k0) = Jz_out(i0, j0+1, k0) + q_w * (dz / dt) * Fz01
    Jz_out(i0+1, j0, k0) = Jz_out(i0+1, j0, k0) + q_w * (dz / dt) * Fz10
    Jz_out(i0+1, j0+1, k0) = Jz_out(i0+1, j0+1, k0) + q_w * (dz / dt) * Fz11
end subroutine Deposit_Segment_VB
!*********************************************************************
subroutine Update_B_FDTD(Bx, By, Bz, Ex, Ey, Ez, dt)      ! 磁场更新器
    use mesh_parameter
    implicit none
    real(8), intent(inout) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:)
    real(8), intent(in) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
    real(8), intent(in) :: dt
    integer :: i, j, k

    do k = 1, nz - 1
        do j = 1, ny - 1
            do i = 1, nx - 1
                Bx(i,j,k) = Bx(i,j,k) - dt*((Ez(i,j+1,k) - Ez(i,j,k)) / dy - (Ey(i,j,k+1) - Ey(i,j,k)) / dz)
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx - 1
                By(i,j,k) = By(i,j,k) - dt*((Ex(i,j,k+1) - Ex(i,j,k)) / dz - (Ez(i+1,j,k) - Ez(i,j,k)) / dx)
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny - 1
            do i = 1, nx - 1
                Bz(i,j,k) = Bz(i,j,k) - dt*((Ey(i+1,j,k) - Ey(i,j,k)) / dx - (Ex(i,j+1,k) - Ex(i,j,k)) / dy)
            end do
        end do
    end do
end subroutine Update_B_FDTD
!*********************************************************************
subroutine Update_E_FDTD(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, dt)      ! 电场更新器
    use mesh_parameter; use physicals_parameter
    implicit none
    real(8), intent(inout) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
    real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:), Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    real(8), intent(in) :: dt
    integer :: i, j, k
    real(8), parameter :: c2 = c*c
    real(8) :: curlB_x, curlB_y, curlB_z

    do k = 1, nz - 1
        do j = 2, ny - 1
            do i = 1, nx - 1
                curlB_x = (Bz(i,j,k) - Bz(i,j-1,k)) / dy - (By(i,j,k) - By(i,j,k-1)) / dz
                Ex(i,j,k) = Ex(i,j,k) + dt * (c2 * curlB_x - Jx(i,j,k)/eps0)
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny - 1
            do i = 1, nx - 1
                curlB_y = (Bx(i,j,k) - Bx(i,j,k-1)) / dz - (Bz(i,j,k) - Bz(i-1,j,k)) / dx
                Ey(i,j,k) = Ey(i,j,k) + dt * (c2 * curlB_y - Jy(i,j,k)/eps0)
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 2, ny - 1
            do i = 1, nx - 1
                curlB_z = (By(i,j,k) - By(i-1,j,k)) / dx - (Bx(i,j,k) - Bx(i,j-1,k)) / dy
                Ez(i,j,k) = Ez(i,j,k) + dt * (c2 * curlB_z - Jz(i,j,k)/eps0)
            end do
        end do
    end do

    ! 边界条件(PEC)
    Ex(:, 1, :) = 0.0d0; Ex(:, ny, :) = 0.0d0; Ez(:, 1, :) = 0.0d0; Ez(:, ny, :) = 0.0d0
end subroutine Update_E_FDTD
!*********************************************************************
subroutine Filter_Currents_3D(Jx, Jy, Jz)   ! 保散度二项式平滑横向平滑
    use mesh_parameter, only: nx, ny, nz
    implicit none
    real(8), intent(inout) :: Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    integer :: pass
    
    do pass = 1, 3
        call Apply_Binomial_Filter(Jx, Jy, Jz)
        call Ghost_Cells(Jx, 'Jx'); call Ghost_Cells(Jy, 'Jy'); call Ghost_Cells(Jz, 'Jz')
    end do
end subroutine Filter_Currents_3D
!*********************************************************************
subroutine Apply_Binomial_Filter(Jx, Jy, Jz)
    use mesh_parameter, only: nx, ny, nz
    implicit none
    real(8), intent(inout) :: Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    real(8), allocatable, save :: temp(:,:,:)
    integer :: i, j, k

    if (.not. allocated(temp)) allocate(temp(0:nx+1, 0:ny+1, 0:nz+1))
    ! Jx在YZ方向横向滤波
    temp = Jx
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                Jx(i,j,k) = 0.5d0 * temp(i,j,k) + 0.25d0 * (temp(i,j-1,k) + temp(i,j+1,k))  !Y平滑
            end do
        end do
    end do
    temp = Jx
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                Jx(i,j,k) = 0.5d0 * temp(i,j,k) + 0.25d0 * (temp(i,j,k-1) + temp(i,j,k+1))  !Z平滑
            end do
        end do
    end do
    ! Jy在XZ方向横向滤波
    temp = Jy
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                Jy(i,j,k) = 0.5d0 * temp(i,j,k) + 0.25d0 * (temp(i-1,j,k) + temp(i+1,j,k))  !X平滑
            end do
        end do
    end do
    temp = Jy
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                Jy(i,j,k) = 0.5d0 * temp(i,j,k) + 0.25d0 * (temp(i,j,k-1) + temp(i,j,k+1))  !Z平滑
            end do
        end do
    end do
    ! Jz在XY方向横向滤波
    temp = Jz
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                Jz(i,j,k) = 0.5d0 * temp(i,j,k) + 0.25d0 * (temp(i-1,j,k) + temp(i+1,j,k))  !X平滑
            end do
        end do
    end do
    temp = Jz
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                Jz(i,j,k) = 0.5d0 * temp(i,j,k) + 0.25d0 * (temp(i,j-1,k) + temp(i,j+1,k))  !Y平滑
            end do
        end do
    end do
end subroutine Apply_Binomial_Filter
!*********************************************************************
subroutine Ghost_Cells(F, var_name)     ! 填充卫星层
    implicit none
    real(8), intent(inout) :: F(0:, 0:, 0:) 
    character(len=*), intent(in) :: var_name
    integer :: nx_local, ny_local, nz_local

    nx_local = ubound(F, 1) - 1; ny_local = ubound(F, 2) - 1; nz_local = ubound(F, 3) - 1

	! X 方向：周期性边界 (连接左右边界)
	F(nx_local, :, :) = F(1, :, :)
	F(0, :, :) = F(nx_local-1, :, :)
	F(nx_local+1, :, :) = F(2, :, :)

    ! Y 方向：理想导电壁边界
    select case (trim(var_name))
        case ('Ex', 'Jx', 'By', 'Bz')
            F(:, 0, :) = -F(:, 1, :)
            F(:, ny_local+1, :) = -F(:, ny_local, :)
        case ('Ey', 'Jy', 'Bx')
            F(:, 0, :) = F(:, 1, :)
            F(:, ny_local+1, :) = F(:, ny_local, :)
        case ('Ez', 'Jz')
            F(:, 0, :) = -F(:, 1, :)
            F(:, ny_local+1, :) = -F(:, ny_local, :)
        case default
            print *, "ERROR: Invalid variable in Ghost_Cells Y-boundary"
            stop
    end select

	! Z 方向：周期性边界 (连接前后边界)
	F(:, :, nz_local) = F(:, :, 1)
	F(:, :, 0) = F(:, :, nz_local-1)
	F(:, :, nz_local+1) = F(:, :, 2)
end subroutine Ghost_Cells
!*********************************************************************
end program pic_equilibrium
