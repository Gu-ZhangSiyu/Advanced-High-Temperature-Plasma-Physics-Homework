!***************************************************************
!   模拟断点续算与重启数据加载模块（位置、速度和二维场）
!
!   调用示例：
!   call load_restart_fields(Bx, By, Bz, Ex, Ey, Ez, restart_step, domain)
!   call read_restart_particles(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, &
!								rank, nprocs, ierr, restart_step, weight_ion, weight_elec)
!   参数设置：
!   (数据结构严格依赖 mesh_parameter 与 pic_types 模块，路径默认为 './for_restart/')
!   
!   传递值说明：
!   Bx,By,Bz,Ex,Ey,Ez		本地重启场 slab，z 向包含两个 ghost 平面
!   ion_local/elec_local	本地进程的离子/电子群组		派生类型粒子数组，恢复位置、速度、ID及权重数据
!   step					重启步长序号				用于定位文件名（如 Bx_1000.txt）的标量整数
!   active_count			探测到的全局活跃粒子总数	用于动态分配内存的整数
!   动态输出不保存逐粒子权重；w_ion/w_elec 为重启时使用的统一权重。
!***************************************************************
module restart_manager_module_CARA
    use mesh_parameter, only: nx, nz, y_plane
    use pic_types, only: particle
    use mpi_f08
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use Domain_Decomposition_Module, only: domain_decomposition
    implicit none

contains
!*********************************************************************
subroutine load_restart_fields(Bx, By, Bz, Ex, Ey, Ez, step, domain)    ! 读取二维文本场并分发本地 slab
	implicit none
	real(8), intent(inout) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:), &
								Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
	integer, intent(in) :: step
	type(domain_decomposition), intent(in) :: domain
	character(len=128) :: filename
	real(8), allocatable :: full_xz(:,:)

	allocate(full_xz(nx, nz))
	write(filename, '("for_restart/Bx_", I0, ".txt")') step
	call read_restart_field_matrix(filename, full_xz)
	Bx(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
	write(filename, '("for_restart/By_", I0, ".txt")') step
	call read_restart_field_matrix(filename, full_xz)
	By(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
	write(filename, '("for_restart/Bz_", I0, ".txt")') step
	call read_restart_field_matrix(filename, full_xz)
	Bz(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
	write(filename, '("for_restart/Ex_", I0, ".txt")') step
	call read_restart_field_matrix(filename, full_xz)
	Ex(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
	write(filename, '("for_restart/Ey_", I0, ".txt")') step
	call read_restart_field_matrix(filename, full_xz)
	Ey(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
	write(filename, '("for_restart/Ez_", I0, ".txt")') step
	call read_restart_field_matrix(filename, full_xz)
	Ez(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
	deallocate(full_xz)
end subroutine load_restart_fields
!*********************************************************************
subroutine read_restart_field_matrix(filename, field)
	character(len=*), intent(in) :: filename
	real(8), intent(out) :: field(:,:)
	integer :: unit_no, io_stat

	open(newunit=unit_no, file=filename, form='formatted', status='old', &
		 action='read', iostat=io_stat)
	if (io_stat /= 0) then
		print *, 'ERROR: cannot open restart field text file: ', trim(filename)
		error stop 33
	end if
	read(unit_no, *, iostat=io_stat) field
	close(unit_no)
	if (io_stat /= 0) then
		print *, 'ERROR: cannot read restart field text matrix: ', trim(filename)
		error stop 34
	end if
end subroutine read_restart_field_matrix
!*********************************************************************
subroutine read_restart_particles(ion_local, elec_local, nIon_local, ion_offset, &
                                  nElec_local, elec_offset, rank, nprocs, ierr, step, w_ion, w_elec)
	implicit none
	type(particle), intent(inout) :: ion_local(:), elec_local(:)
	integer, intent(in) :: nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, step
	integer, intent(out) :: ierr
	real(8), intent(in) :: w_ion, w_elec
	call read_restart_species(ion_local, nIon_local, ion_offset, rank, nprocs, ierr, 'ions', step, w_ion)
	call read_restart_species(elec_local, nElec_local, elec_offset, rank, nprocs, ierr, 'electrons', step, w_elec)
end subroutine read_restart_particles
!*********************************************************************	
subroutine read_restart_species(part, nLocal, offset, rank, nprocs, ierr, species_in, step, macro_weight)	! 重启读取粒子 (MPI-IO)
	implicit none
	type(particle), intent(inout) :: part(:)
	integer, intent(in) :: nLocal, offset, rank, nprocs, step
	integer, intent(out) :: ierr
	character(len=*), intent(in) :: species_in
	real(8), intent(in) :: macro_weight
	character(len=64) :: fname_pos, fname_vel
	character(len=12) :: species
	real(8), allocatable :: data_pos(:,:), data_vel(:,:)
	integer :: i
    integer(kind=8) :: nGlobal
	integer(kind=MPI_OFFSET_KIND) :: disp
	type(MPI_File) :: fh_pos, fh_vel
	type(MPI_Status) :: status

	species = adjustl(species_in)
	write(fname_pos, '("for_restart/",A,"_position_",I0,".bin")') trim(species), step
	write(fname_vel, '("for_restart/",A,"_velocity_",I0,".bin")') trim(species), step
    call get_active_particle_count(species_in, step, nGlobal)
    call require_restart_file(fname_vel, nGlobal * 32_8)
	if (.not. ieee_is_finite(macro_weight) .or. macro_weight <= 0.0d0) then
        print *, 'ERROR: restart uniform weight must be positive for ', trim(species)
        call MPI_Abort(MPI_COMM_WORLD, 44, ierr)
    end if
    allocate(data_pos(4, max(nLocal, 1)), data_vel(4, max(nLocal, 1)))
	disp = int(offset, MPI_OFFSET_KIND) * 32_MPI_OFFSET_KIND
	call MPI_File_open(MPI_COMM_WORLD, trim(fname_pos), MPI_MODE_RDONLY, MPI_INFO_NULL, fh_pos, ierr)
	if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 45, ierr)
	call MPI_File_read_at_all(fh_pos, disp, data_pos(1,1), nLocal * 4, MPI_DOUBLE_PRECISION, status, ierr)
	call check_restart_read(status, ierr, nLocal*4)
	call MPI_File_close(fh_pos, ierr)
	call MPI_File_open(MPI_COMM_WORLD, trim(fname_vel), MPI_MODE_RDONLY, MPI_INFO_NULL, fh_vel, ierr)
	if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 45, ierr)
	call MPI_File_read_at_all(fh_vel, disp, data_vel(1,1), nLocal * 4, MPI_DOUBLE_PRECISION, status, ierr)
    call check_restart_read(status, ierr, nLocal*4)
	call MPI_File_close(fh_vel, ierr)
	    if (.not. all(ieee_is_finite(data_pos(:,1:nLocal))) .or. &
	        .not. all(ieee_is_finite(data_vel(:,1:nLocal)))) call MPI_Abort(MPI_COMM_WORLD, 46, ierr)
    if (any(data_pos(1,1:nLocal) /= data_vel(1,1:nLocal))) call MPI_Abort(MPI_COMM_WORLD, 46, ierr)
	    if (any(data_pos(1,1:nLocal) < 1.0d0) .or. &
	        any(data_pos(1,1:nLocal) > 9007199254740991.0d0)) call MPI_Abort(MPI_COMM_WORLD, 46, ierr)
	do i = 1, nLocal
		part(i)%id = int(data_pos(1, i), 8)
        if (dble(part(i)%id) /= data_pos(1,i)) call MPI_Abort(MPI_COMM_WORLD, 46, ierr)
		part(i)%x = data_pos(2, i)
		part(i)%y = y_plane
		part(i)%z = data_pos(4, i)
		part(i)%vx = data_vel(2, i)
		part(i)%vy = data_vel(3, i)
		part(i)%vz = data_vel(4, i)
		part(i)%weight = macro_weight
	end do
	deallocate(data_pos, data_vel)
end subroutine read_restart_species
!*********************************************************************
subroutine get_active_particle_count(species_in, step, active_count)
    character(len=*), intent(in) :: species_in
    integer, intent(in) :: step
    integer(kind=8), intent(out) :: active_count
    character(len=128) :: filename
    integer(kind=8) :: file_bytes
    integer :: ios, ierr
    logical :: exists

    write(filename, '("for_restart/",A,"_position_",I0,".bin")') trim(species_in), step
    inquire(file=filename, exist=exists, size=file_bytes, iostat=ios)
    if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 47, ierr)
    if (.not. exists) call MPI_Abort(MPI_COMM_WORLD, 47, ierr)
    if (file_bytes < 0_8) call MPI_Abort(MPI_COMM_WORLD, 47, ierr)
    if (mod(file_bytes, 32_8) /= 0_8) call MPI_Abort(MPI_COMM_WORLD, 47, ierr)
    active_count = file_bytes / 32_8
end subroutine get_active_particle_count
!*********************************************************************
subroutine require_restart_file(filename, expected_bytes)
    character(len=*), intent(in) :: filename
    integer(kind=8), intent(in) :: expected_bytes
    integer(kind=8) :: file_bytes
    integer :: ios, ierr
    logical :: exists
    inquire(file=filename, exist=exists, size=file_bytes, iostat=ios)
    if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 48, ierr)
    if (.not. exists) call MPI_Abort(MPI_COMM_WORLD, 48, ierr)
    if (file_bytes /= expected_bytes) call MPI_Abort(MPI_COMM_WORLD, 48, ierr)
end subroutine require_restart_file
!*********************************************************************
subroutine check_restart_read(status, read_error, expected)
    type(MPI_Status), intent(in) :: status
    integer, intent(in) :: read_error, expected
    integer :: actual, ierr
    if (read_error /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 49, ierr)
    call MPI_Get_count(status, MPI_DOUBLE_PRECISION, actual, ierr)
    if (actual /= expected) call MPI_Abort(MPI_COMM_WORLD, 49, ierr)
end subroutine check_restart_read
!*********************************************************************
subroutine initialize_injection_state(part, nLocal, last_id, fraction)
    implicit none
    type(particle), intent(in) :: part(:)
    integer, intent(in) :: nLocal
    integer(kind=8), intent(out) :: last_id
    real(8), intent(out) :: fraction
    integer(kind=8) :: local_max
    integer :: ierr

    ! 动态state文件已取消；重启时以当前检查点中的最大ID作为新的注入起点。
    local_max = 0_8
    if (nLocal > 0) local_max = maxval(part(1:nLocal)%id)
    call MPI_Allreduce(local_max, last_id, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr)
    fraction = 0.0d0
end subroutine initialize_injection_state
!*********************************************************************
end module restart_manager_module_CARA
