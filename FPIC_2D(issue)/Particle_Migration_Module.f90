!*********************************************************************
! 粒子空间归属、动态缓冲区和 MPI 迁移
!*********************************************************************
module Particle_Migration_Module
    use mpi_f08
    use pic_types, only: particle
    use Domain_Decomposition_Module, only: domain_decomposition, Owner_Of_Position
    implicit none
    private
    integer, parameter :: packed_particle_size = 8
    public :: Ensure_Particle_Capacity, Redistribute_Particles_Global
    public :: Pack_Particle, Unpack_Particle, packed_particle_size
    public :: Migrate_Particles_Neighbors, Remove_Outside_Particles, Report_Particle_Load

contains
!*********************************************************************
subroutine Ensure_Particle_Capacity(part, n_active, capacity, required)
    type(particle), pointer, intent(inout) :: part(:)
    integer, intent(in) :: n_active, required
    integer, intent(inout) :: capacity
    type(particle), pointer :: grown(:)
    integer :: new_capacity
    if (required <= capacity) return
    new_capacity = max(required, max(1024, capacity + max(capacity / 2, 1024)))
    allocate(grown(new_capacity))
    if (n_active > 0) grown(1:n_active) = part(1:n_active)
    if (associated(part)) deallocate(part)
    part => grown; capacity = new_capacity
end subroutine Ensure_Particle_Capacity
!*********************************************************************
subroutine Pack_Particle(p, buffer)
    type(particle), intent(in) :: p
    real(8), intent(out) :: buffer(packed_particle_size)
    buffer(1) = transfer(p%id, buffer(1))
    buffer(2) = p%x; buffer(3) = p%y; buffer(4) = p%z
    buffer(5) = p%vx; buffer(6) = p%vy; buffer(7) = p%vz
    buffer(8) = p%weight
end subroutine Pack_Particle
!*********************************************************************
subroutine Unpack_Particle(buffer, p)
    real(8), intent(in) :: buffer(packed_particle_size)
    type(particle), intent(out) :: p
    p%id = transfer(buffer(1), p%id)
    p%x = buffer(2); p%y = buffer(3); p%z = buffer(4)
    p%vx = buffer(5); p%vy = buffer(6); p%vz = buffer(7)
    p%weight = buffer(8)
end subroutine Unpack_Particle
!*********************************************************************
subroutine Redistribute_Particles_Global(part, n_active, capacity, domain)
    type(particle), pointer, intent(inout) :: part(:)
    integer, intent(inout) :: n_active, capacity
    type(domain_decomposition), intent(in) :: domain
    integer, allocatable :: send_particles(:), recv_particles(:), cursor(:)
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
    real(8), allocatable :: send_buffer(:), recv_buffer(:)
    integer :: i, dest, total_send, total_recv, ierr, slot
    allocate(send_particles(0:domain%nprocs-1), recv_particles(0:domain%nprocs-1))
    allocate(send_counts(0:domain%nprocs-1), recv_counts(0:domain%nprocs-1))
    allocate(send_displs(0:domain%nprocs-1), recv_displs(0:domain%nprocs-1))
    allocate(cursor(0:domain%nprocs-1))

    send_particles = 0
    do i = 1, n_active
        if (part(i)%weight <= 0.0d0) cycle
        dest = Owner_Of_Position(part(i)%z, domain)
        if (dest /= MPI_PROC_NULL) send_particles(dest) = send_particles(dest) + 1
    end do
    call MPI_Alltoall(send_particles, 1, MPI_INTEGER, recv_particles, 1, MPI_INTEGER, domain%comm, ierr)
    send_counts = packed_particle_size * send_particles; recv_counts = packed_particle_size * recv_particles
    send_displs(0) = 0; recv_displs(0) = 0
    do i = 1, domain%nprocs - 1
        send_displs(i) = send_displs(i-1) + send_counts(i-1)
        recv_displs(i) = recv_displs(i-1) + recv_counts(i-1)
    end do
    total_send = sum(send_particles); total_recv = sum(recv_particles)
    allocate(send_buffer(max(1, packed_particle_size * total_send)))
    allocate(recv_buffer(max(1, packed_particle_size * total_recv)))
    cursor = send_displs
	
    do i = 1, n_active
        if (part(i)%weight <= 0.0d0) cycle
        dest = Owner_Of_Position(part(i)%z, domain)
        if (dest == MPI_PROC_NULL) cycle
        slot = cursor(dest) + 1
        call Pack_Particle(part(i), send_buffer(slot:slot+packed_particle_size-1))
        cursor(dest) = cursor(dest) + packed_particle_size
    end do
    call MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_DOUBLE_PRECISION, &
                       recv_buffer, recv_counts, recv_displs, MPI_DOUBLE_PRECISION, domain%comm, ierr)
    call Ensure_Particle_Capacity(part, 0, capacity, total_recv)
    do i = 1, total_recv
        slot = (i - 1) * packed_particle_size + 1
        call Unpack_Particle(recv_buffer(slot:slot+packed_particle_size-1), part(i))
    end do
    n_active = total_recv

    deallocate(send_particles, recv_particles, send_counts, recv_counts)
    deallocate(send_displs, recv_displs, cursor, send_buffer, recv_buffer)
end subroutine Redistribute_Particles_Global
!*********************************************************************
subroutine Migrate_Particles_Neighbors(part, n_active, capacity, domain, species_name)
    type(particle), pointer, intent(inout) :: part(:)
    integer, intent(inout) :: n_active, capacity
    type(domain_decomposition), intent(in) :: domain
    character(len=*), intent(in) :: species_name
    real(8), allocatable :: send_left(:,:), send_right(:,:), recv_left(:,:), recv_right(:,:)
    type(MPI_Status) :: status
    integer :: i, owner, keep_count, n_send_left, n_send_right
    integer :: n_recv_left, n_recv_right, ierr, required

    n_send_left = 0; n_send_right = 0; keep_count = 0
    do i = 1, n_active		! 第一遍只计数，避免每一步按全部粒子数申请两个发送缓冲区
        if (part(i)%weight <= 0.0d0) cycle
        owner = Owner_Of_Position(part(i)%z, domain)
        if (owner == domain%rank .or. owner == MPI_PROC_NULL) then
            keep_count = keep_count + 1
        else if (owner == domain%left) then
            n_send_left = n_send_left + 1
        else if (owner == domain%right) then
            n_send_right = n_send_right + 1
        else
            print *, 'ERROR: ', trim(species_name), ' particle crossed a non-neighbor domain on rank ', domain%rank
            print *, 'particle id=', part(i)%id, ' z=', part(i)%z, ' owner=', owner
            call MPI_Abort(domain%comm, 21, ierr)
        end if
    end do

    allocate(send_left(packed_particle_size, max(1, n_send_left)))
    allocate(send_right(packed_particle_size, max(1, n_send_right)))
    n_send_left = 0; n_send_right = 0; keep_count = 0
    do i = 1, n_active
        if (part(i)%weight <= 0.0d0) cycle
        owner = Owner_Of_Position(part(i)%z, domain)
        if (owner == domain%rank .or. owner == MPI_PROC_NULL) then
            keep_count = keep_count + 1
            if (keep_count /= i) part(keep_count) = part(i)
        else if (owner == domain%left) then
            n_send_left = n_send_left + 1
            call Pack_Particle(part(i), send_left(:,n_send_left))
        else if (owner == domain%right) then
            n_send_right = n_send_right + 1
            call Pack_Particle(part(i), send_right(:,n_send_right))
        end if
    end do

    n_recv_left = 0; n_recv_right = 0
    call MPI_Sendrecv(n_send_left, 1, MPI_INTEGER, domain%left, 801, &
                      n_recv_right, 1, MPI_INTEGER, domain%right, 801, domain%comm, status, ierr)
    call MPI_Sendrecv(n_send_right, 1, MPI_INTEGER, domain%right, 802, &
                      n_recv_left, 1, MPI_INTEGER, domain%left, 802, domain%comm, status, ierr)
    allocate(recv_left(packed_particle_size, max(1, n_recv_left)))
    allocate(recv_right(packed_particle_size, max(1, n_recv_right)))
    call MPI_Sendrecv(send_left, packed_particle_size*n_send_left, MPI_DOUBLE_PRECISION, domain%left, 803, &
                      recv_right, packed_particle_size*n_recv_right, MPI_DOUBLE_PRECISION, domain%right, 803, domain%comm, status, ierr)
    call MPI_Sendrecv(send_right, packed_particle_size*n_send_right, MPI_DOUBLE_PRECISION, domain%right, 804, &
                      recv_left, packed_particle_size*n_recv_left, MPI_DOUBLE_PRECISION, domain%left, 804, domain%comm, status, ierr)

    required = keep_count + n_recv_left + n_recv_right
    call Ensure_Particle_Capacity(part, keep_count, capacity, required)
    do i = 1, n_recv_left; call Unpack_Particle(recv_left(:,i), part(keep_count+i)); end do
    do i = 1, n_recv_right; call Unpack_Particle(recv_right(:,i), part(keep_count+n_recv_left+i)); end do
	n_active = required

    deallocate(send_left, send_right, recv_left, recv_right)
end subroutine Migrate_Particles_Neighbors
!*********************************************************************
subroutine Remove_Outside_Particles(part, n_active)
    use mesh_parameter, only: x_min, x_max, z_min, z_max, y_plane
    type(particle), intent(inout) :: part(:)
    integer, intent(inout) :: n_active
    integer :: i, keep_count

    keep_count = 0
    do i = 1, n_active
        if (part(i)%weight <= 0.0d0) cycle
        if (part(i)%x < x_min .or. part(i)%x >= x_max) cycle
        if (part(i)%z < z_min .or. part(i)%z >= z_max) cycle
        part(i)%y = y_plane
        keep_count = keep_count + 1
        if (keep_count /= i) part(keep_count) = part(i)
    end do
    n_active = keep_count
end subroutine Remove_Outside_Particles

subroutine Report_Particle_Load(n_ion, n_elec, step, domain)
    integer, intent(in) :: n_ion, n_elec, step
    type(domain_decomposition), intent(in) :: domain
    integer :: local_count, min_count, max_count, sum_count, ierr
    real(8) :: imbalance

    local_count = n_ion + n_elec
    call MPI_Reduce(local_count, min_count, 1, MPI_INTEGER, MPI_MIN, 0, domain%comm, ierr)
    call MPI_Reduce(local_count, max_count, 1, MPI_INTEGER, MPI_MAX, 0, domain%comm, ierr)
    call MPI_Reduce(local_count, sum_count, 1, MPI_INTEGER, MPI_SUM, 0, domain%comm, ierr)
    if (domain%rank == 0) then
        if (sum_count > 0) then
            imbalance = dble(max_count) / (dble(sum_count) / dble(domain%nprocs))
        else
            imbalance = 1.0d0
        end if
        print '(A,I10,A,I12,A,I12,A,I12,A,F8.3)', 'load step=', step, &
              ' min=', min_count, ' max=', max_count, ' total=', sum_count, &
              ' max/avg=', imbalance
    end if
end subroutine Report_Particle_Load

end module Particle_Migration_Module
