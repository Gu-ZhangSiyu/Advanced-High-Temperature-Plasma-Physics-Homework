!***************************************************************
!   PIC背景平衡场与粒子数据并行流式读取子程序
!
!   调用示例
!   call load_equilibrium(Bx0, By0, Bz0, Ex0, Ey0, Ez0, domain)
!   call read_particle_data(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, ierr)
!
!   参数设置：
!   (无直接硬编码控制参数，网格尺寸与基础数据结构严格依赖 mesh_parameter 与 pic_types 模块)
!   
!   传递值说明
!   Bx,By,Bz,Ex,Ey,Ez			本地背景场 slab，z 向包含两个 ghost 平面
!
!   传递值说明
!   ion_local/elec_local		本地进程的离子/电子群组			派生类型粒子数组，装载相空间坐标与权重数据
!   nIon_local/nElec_local		本地进程分配的离子/电子数		标量值，决定当前MPI rank的本地内存容量
!   ion_offset/elec_offset		全局粒子ID的分布偏移量			标量值，计算底层二进制文件寻址(byte_offset)的核心基准
!   rank						当前MPI进程编号				标量值，用于并行分块I/O的进程识别
!   nprocs						MPI通讯域内总进程数			标量值，用于 Allgather 规约全局拓扑
!   ierr						MPI标准错误返回码				标量值
!***************************************************************
module data_loader_module
    use mesh_parameter, only: nx, nz, y_plane
    use pic_types, only: particle
    use mpi_f08
    use Domain_Decomposition_Module, only: domain_decomposition
    implicit none
    private
    public :: load_equilibrium, read_particle_data

contains
!*********************************************************************
subroutine load_equilibrium(Bx, By, Bz, Ex, Ey, Ez, domain)		! 读取二维 XZ 背景场文本矩阵并分发本地 slab
    implicit none
    real(8), intent(out) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:), Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    real(8), allocatable :: full_xz(:,:)

    allocate(full_xz(nx, nz))
    call read_field_matrix_text('out_xz/Bx_xz.txt', full_xz)
    Bx(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
    call read_field_matrix_text('out_xz/By_xz.txt', full_xz)
    By(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
    call read_field_matrix_text('out_xz/Bz_xz.txt', full_xz)
    Bz(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
    call read_field_matrix_text('out_xz/Ex_xz.txt', full_xz)
    Ex(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
    call read_field_matrix_text('out_xz/Ey_xz.txt', full_xz)
    Ey(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
    call read_field_matrix_text('out_xz/Ez_xz.txt', full_xz)
    Ez(1:nx,1,1:domain%nz_local) = full_xz(:,domain%k_start:domain%k_start+domain%nz_local-1)
    deallocate(full_xz)
end subroutine load_equilibrium
!*********************************************************************
subroutine read_field_matrix_text(filename, field)
    character(len=*), intent(in) :: filename
    real(8), intent(out) :: field(:,:)
    integer :: unit_no, ios, ierr
    logical :: exists

    inquire(file=filename, exist=exists, iostat=ios)
    if (ios /= 0 .or. .not. exists) then
        print *, 'ERROR: cannot find XZ field text file: ', trim(filename)
        call MPI_Abort(MPI_COMM_WORLD, 36, ierr)
    end if
    open(newunit=unit_no, file=filename, form='formatted', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'ERROR: cannot open XZ field text file: ', trim(filename)
        call MPI_Abort(MPI_COMM_WORLD, 37, ierr)
    end if
    read(unit_no, *, iostat=ios) field
    close(unit_no)
    if (ios /= 0) then
        print *, 'ERROR: cannot read XZ field text matrix: ', trim(filename)
        call MPI_Abort(MPI_COMM_WORLD, 38, ierr)
    end if
end subroutine read_field_matrix_text
!*********************************************************************
subroutine read_particle_data(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, ierr)
    use control_parameter, only: nIon, nElec, background_weight_ion, background_weight_elec
    type(particle), intent(out) :: ion_local(:), elec_local(:)
    integer, intent(in) :: nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs
    integer, intent(out) :: ierr
    call read_background_species(ion_local, nIon_local, ion_offset, nIon, 'ion', background_weight_ion)
    call read_background_species(elec_local, nElec_local, elec_offset, nElec, 'elec', background_weight_elec)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, 'Background particle loading complete.'
end subroutine read_particle_data
!*********************************************************************
subroutine read_background_species(part, nLocal, offset, nGlobal, species, manual_weight)
    use control_parameter, only: background_weights_from_file
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(particle), intent(out) :: part(:)
    integer, intent(in) :: nLocal, offset
    integer(kind=8), intent(in) :: nGlobal
    character(len=*), intent(in) :: species
    real(8), intent(in) :: manual_weight
    ! 2D3V 背景粒子只读取 x,z 位置；y 位置固定为 y_plane。
    character(len=6), parameter :: components(6) = [character(len=6) :: 'x', 'z', 'vx', 'vy', 'vz', 'weight']
    character(len=128) :: filename
    real(8), allocatable :: values(:)
    integer :: c, nc, i, u, ios, ierr
    integer(kind=8) :: file_bytes
    logical :: exists

    if (nGlobal == 0_8) return
    nc = 5
    if (background_weights_from_file) nc = 6
    allocate(values(nLocal))
    if (.not. background_weights_from_file) then
        if (.not. ieee_is_finite(manual_weight) .or. manual_weight <= 0.0d0) then
            print *, 'ERROR: positive background weight required for ', species
            call MPI_Abort(MPI_COMM_WORLD, 35, ierr)
        end if
        if (nLocal > 0) part(1:nLocal)%weight = manual_weight
    end if
    if (nLocal > 0) part(1:nLocal)%y = y_plane
    do c = 1, nc
        filename = 'init_particles/' // species // '_' // trim(components(c)) // '.bin'
        inquire(file=filename, exist=exists, size=file_bytes, iostat=ios)
        if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 36, ierr)
        if (.not. exists) call MPI_Abort(MPI_COMM_WORLD, 36, ierr)
        if (file_bytes /= nGlobal * 8_8) then
            print *, 'ERROR: background file count mismatch: ', trim(filename)
            call MPI_Abort(MPI_COMM_WORLD, 36, ierr)
        end if
        if (nLocal == 0) cycle
        open(newunit=u, file=filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 37, ierr)
        read(u, pos=int(offset,8)*8_8+1_8, iostat=ios) values
        close(u)
        if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 38, ierr)
        if (.not. all(ieee_is_finite(values))) then
            print *, 'ERROR: non-finite background data: ', trim(filename)
            call MPI_Abort(MPI_COMM_WORLD, 39, ierr)
        end if
        select case (c)
        case (1); part(1:nLocal)%x = values; case (2); part(1:nLocal)%z = values
        case (3); part(1:nLocal)%vx = values; case (4); part(1:nLocal)%vy = values
        case (5); part(1:nLocal)%vz = values
        case (6)
            if (any(values < 0.0d0)) call MPI_Abort(MPI_COMM_WORLD, 39, ierr)
            part(1:nLocal)%weight = values
        end select
    end do
    do i = 1, nLocal; part(i)%id = int(offset,8) + int(i,8); end do
    deallocate(values)
end subroutine read_background_species
!*********************************************************************
end module data_loader_module
