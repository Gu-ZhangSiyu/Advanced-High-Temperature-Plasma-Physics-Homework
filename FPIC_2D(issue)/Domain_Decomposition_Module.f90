!*********************************************************************
!	一维z向区域分解, 场halo交换与电流边界累加
!*********************************************************************
module Domain_Decomposition_Module
    use mpi_f08
    use mesh_parameter, only: nx, ny, nz, z_min, z_max, dz, inv_dz
    implicit none
    private
!*********************************************************************
    type, public :: domain_decomposition
        type(MPI_Comm) :: comm
        integer :: rank = -1
        integer :: nprocs = 0
        integer, allocatable :: cuts(:) ! 0:nprocs, global starts; last entry nz+1
        integer :: nz_local = 0
        integer :: k_start = 0		! 本rank拥有的第一个全局z节点(1-based)
        integer :: k_end = -1		! 本rank拥有的最后一个全局z节点
        integer :: left = MPI_PROC_NULL
        integer :: right = MPI_PROC_NULL
        real(8) :: z_origin = 0.0d0	! 本地数组k=1对应的物理坐标
    end type domain_decomposition
    public :: Initialize_Domain, Owner_Of_Position, Exchange_Field_Halos
    public :: Accumulate_Current_Halos, Report_Domain_Decomposition
    public :: Set_Domain_Cuts, Remap_Field

contains
!*********************************************************************
subroutine Initialize_Domain(domain, comm)
    type(domain_decomposition), intent(out) :: domain
    type(MPI_Comm), intent(in) :: comm
    integer :: base, remainder, ierr, min_local, r

    domain%comm = comm
    call MPI_Comm_rank(comm, domain%rank, ierr)
    call MPI_Comm_size(comm, domain%nprocs, ierr)
    if (domain%nprocs > nz / 2) then
        if (domain%rank == 0) then
            print *, 'ERROR: z decomposition requires at least two owned z nodes per rank.'
            print *, '       Reduce MPI ranks to <= ', nz / 2
        end if
        call MPI_Abort(comm, 11, ierr)
    end if
    allocate(domain%cuts(0:domain%nprocs))
    do r = 0, domain%nprocs
        domain%cuts(r) = 1 + r * (nz / domain%nprocs) + min(r, mod(nz, domain%nprocs))
    end do
    base = nz / domain%nprocs
    remainder = mod(nz, domain%nprocs)
    domain%nz_local = base
    if (domain%rank < remainder) domain%nz_local = domain%nz_local + 1
    if (domain%rank < remainder) then
        domain%k_start = domain%rank * (base + 1) + 1
    else
        domain%k_start = remainder * (base + 1) + (domain%rank - remainder) * base + 1
    end if
    domain%k_end = domain%k_start + domain%nz_local - 1
    domain%z_origin = z_min + dble(domain%k_start - 1) * dz
    if (domain%rank > 0) domain%left = domain%rank - 1
    if (domain%rank < domain%nprocs - 1) domain%right = domain%rank + 1
    call MPI_Allreduce(domain%nz_local, min_local, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (min_local < 2) then
        if (domain%rank == 0) print *, 'ERROR: invalid z decomposition (local nz < 2).'
        call MPI_Abort(comm, 12, ierr)
    end if
end subroutine Initialize_Domain
!*********************************************************************
integer function Owner_Of_Position(z, domain) result(owner)
    real(8), intent(in) :: z
    type(domain_decomposition), intent(in) :: domain
    integer :: global_node, lo, hi, mid
    owner = MPI_PROC_NULL
    if (z < z_min .or. z >= z_max) return
    global_node = max(1, min(floor((z - z_min) * inv_dz) + 1, nz - 1))
    lo = 0; hi = domain%nprocs
    do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (global_node < domain%cuts(mid)) then
            hi = mid
        else
            lo = mid
        end if
    end do
    owner = lo
end function Owner_Of_Position
!*********************************************************************
! 从同一个全局切割表中设置rank的所有权
subroutine Set_Domain_Cuts(domain, cuts)
    type(domain_decomposition), intent(inout) :: domain
    integer, intent(in) :: cuts(0:)
    integer :: ierr
    if (size(cuts) /= domain%nprocs + 1) call MPI_Abort(domain%comm, 13, ierr)
    if (cuts(0) /= 1 .or. cuts(domain%nprocs) /= nz+1) call MPI_Abort(domain%comm, 13, ierr)
    if (any(cuts(1:) - cuts(:domain%nprocs-1) < 2)) call MPI_Abort(domain%comm, 13, ierr)
    domain%cuts = cuts
    domain%k_start = cuts(domain%rank)
    domain%k_end = cuts(domain%rank+1) - 1
    domain%nz_local = domain%k_end - domain%k_start + 1
    domain%z_origin = z_min + dble(domain%k_start - 1) * dz
end subroutine Set_Domain_Cuts
!*********************************************************************
! 将拥有的planes直接在原所有者和新所有者之间转移, 调用方重建特定组件的物理/通信halos
subroutine Remap_Field(F, old_domain, new_domain)
    real(8), allocatable, intent(inout) :: F(:,:,:)
    type(domain_decomposition), intent(in) :: old_domain, new_domain
    real(8), allocatable :: moved(:,:,:)
    integer :: sc(0:old_domain%nprocs-1), rc(0:old_domain%nprocs-1)
    integer :: sd(0:old_domain%nprocs-1), rd(0:old_domain%nprocs-1)
    integer :: r, lo, hi, plane_size, ierr
    plane_size = size(F,1) * size(F,2)
    allocate(moved(0:size(F,1)-1, 0:size(F,2)-1, 0:new_domain%nz_local+1))
    moved = 0.0d0
    sc=0; rc=0; sd=0; rd=0
    do r=0, old_domain%nprocs-1
        lo=max(old_domain%k_start,new_domain%cuts(r)); hi=min(old_domain%k_end,new_domain%cuts(r+1)-1)
        if (lo <= hi) then
            sc(r)=plane_size*(hi-lo+1); sd(r)=plane_size*(lo-old_domain%k_start)
        end if
        lo=max(new_domain%k_start,old_domain%cuts(r)); hi=min(new_domain%k_end,old_domain%cuts(r+1)-1)
        if (lo <= hi) then
            rc(r)=plane_size*(hi-lo+1); rd(r)=plane_size*(lo-new_domain%k_start)
        end if
    end do
    call MPI_Alltoallv(F(:,:,1:old_domain%nz_local), sc, sd, MPI_DOUBLE_PRECISION, &
                      moved(:,:,1:new_domain%nz_local), rc, rd, MPI_DOUBLE_PRECISION, new_domain%comm, ierr)
    call move_alloc(moved,F)
end subroutine Remap_Field
!*********************************************************************
subroutine Exchange_Field_Halos(F, domain)
    real(8), intent(inout), contiguous :: F(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    type(MPI_Status) :: status
    integer :: ierr, plane_size

    plane_size = size(F, 1) * size(F, 2)

    ! 向左发送第一个owned平面,从右邻居接收其第一个owned平面
    call MPI_Sendrecv(F(:,:,1), plane_size, MPI_DOUBLE_PRECISION, domain%left, 701, &
                      F(:,:,domain%nz_local+1), plane_size, MPI_DOUBLE_PRECISION, domain%right, 701, domain%comm, status, ierr)
    ! 向右发送最后一个owned平面,从左邻居接收其最后一个owned平面
    call MPI_Sendrecv(F(:,:,domain%nz_local), plane_size, MPI_DOUBLE_PRECISION, domain%right, 702, &
                      F(:,:,0), plane_size, MPI_DOUBLE_PRECISION, domain%left, 702, domain%comm, status, ierr)
    ! 仅全局物理端点使用原代码的零梯度开放边界
    if (domain%left == MPI_PROC_NULL) F(:,:,0) = F(:,:,1)
    if (domain%right == MPI_PROC_NULL) F(:,:,domain%nz_local+1) = F(:,:,domain%nz_local)
end subroutine Exchange_Field_Halos
!*********************************************************************
subroutine Accumulate_Current_Halos(F, domain)
    real(8), intent(inout), contiguous :: F(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    type(MPI_Status) :: status
    real(8), allocatable :: recv_left(:,:), recv_right(:,:)
    integer :: ierr, plane_size

    plane_size = size(F, 1) * size(F, 2)
    allocate(recv_left(lbound(F,1):ubound(F,1), lbound(F,2):ubound(F,2)))
    allocate(recv_right(lbound(F,1):ubound(F,1), lbound(F,2):ubound(F,2)))
    recv_left = 0.0d0; recv_right = 0.0d0

    ! 下侧ghost上的沉积属于左邻居,同时接收右邻居下侧ghost的贡献
    call MPI_Sendrecv(F(:,:,0), plane_size, MPI_DOUBLE_PRECISION, domain%left, 711, &
                      recv_right, plane_size, MPI_DOUBLE_PRECISION, domain%right, 711, domain%comm, status, ierr)
    ! 上侧ghost上的沉积属于右邻居,同时接收左邻居上侧ghost的贡献
    call MPI_Sendrecv(F(:,:,domain%nz_local+1), plane_size, MPI_DOUBLE_PRECISION, domain%right, 712, &
                      recv_left, plane_size, MPI_DOUBLE_PRECISION, domain%left, 712, domain%comm, status, ierr)
    if (domain%left /= MPI_PROC_NULL) F(:,:,1) = F(:,:,1) + recv_left
    if (domain%right /= MPI_PROC_NULL) F(:,:,domain%nz_local) = F(:,:,domain%nz_local) + recv_right
    F(:,:,0) = 0.0d0; F(:,:,domain%nz_local+1) = 0.0d0
    deallocate(recv_left, recv_right)
end subroutine Accumulate_Current_Halos
!*********************************************************************
subroutine Report_Domain_Decomposition(domain)
    type(domain_decomposition), intent(in) :: domain
    integer :: ierr
    if (domain%rank == 0) then
        print *, 'z-slab domain decomposition enabled: ', domain%nprocs, ' MPI ranks'
        print *, 'global grid = ', nx, ' x ', ny, ' x ', nz
    end if
    call MPI_Barrier(domain%comm, ierr)
    print '(A,I5,A,I6,A,I6,A,I6)', 'rank ', domain%rank, ': global k=', &
          domain%k_start, ' ... ', domain%k_end, ', local nz=', domain%nz_local
    call MPI_Barrier(domain%comm, ierr)
end subroutine Report_Domain_Decomposition
!*********************************************************************
end module Domain_Decomposition_Module
