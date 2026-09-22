! 2D XZ Harris 背景粒子初始化器
!
! 原始程序在 (257,257,3) 网格上按 XY 平面采样。当前工程是 XZ 平面，
! 因此本程序使用 (nx,ny,nz)=(129,1,129)，固定 y=0，保留三分量速度。
! 输出文件名严格匹配 data_loader_module：
!   init_particles/{ion,elec}_{x,z,vx,vy,vz,weight}.bin
!
module mesh_parameter
    implicit none
    integer, parameter :: nx = 129, ny = 1, nz = 129
    real(8), parameter :: x_min = -2.0d-2, x_max =  2.0d-2
    real(8), parameter :: y_min = -2.0d-2, y_max =  2.0d-2
    real(8), parameter :: y_plane = 0.0d0
    real(8), parameter :: z_min = -2.0d-2, z_max =  2.0d-2
    real(8), parameter :: dx = (x_max - x_min) / dble(nx - 1)
    ! ny=1 时 dy 不再是网格间距，而是二维模型的有效出平面厚度。
    real(8), parameter :: dy = y_max - y_min
    real(8), parameter :: dz = (z_max - z_min) / dble(nz - 1)
end module mesh_parameter
!*********************************************************************
module physicals_parameter
    implicit none
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: q = 1.60217733d-19
end module physicals_parameter
!*********************************************************************
module pic_types
    use physicals_parameter, only: q, pi
    implicit none
    type :: particle
        real(8) :: x, y, z
        real(8) :: vx, vy, vz
        real(8) :: charge, mass
        real(8) :: weight
    end type particle
contains
    function sample_maxwell(temp_ev, mass) result(v)
        implicit none
        real(8), intent(in) :: temp_ev, mass
        real(8) :: v(3), sigma, u1, u2, u3, u4

        sigma = sqrt(temp_ev * q / mass)
        call random_number(u1); call random_number(u2)
        call random_number(u3); call random_number(u4)
        u1 = max(u1, 1.0d-14)
        u3 = max(u3, 1.0d-14)
        v(1) = sigma * sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
        v(2) = sigma * sqrt(-2.0d0*log(u1)) * sin(2.0d0*pi*u2)
        v(3) = sigma * sqrt(-2.0d0*log(u3)) * cos(2.0d0*pi*u4)
    end function sample_maxwell
end module pic_types
!*********************************************************************
program pic_equilibrium
    use mpi
    use mesh_parameter
    use physicals_parameter
    use pic_types
    implicit none

    integer :: rank, nprocs, ierr
    integer :: nIon = 10000000
    integer :: nElec = 10000000
    integer :: nIon_local, nElec_local
    integer :: ion_offset, elec_offset
    integer :: base_ion, extra_ion, base_elec, extra_elec
    integer, parameter :: base_seed = 13579
    type(particle), allocatable :: ion_local(:), elec_local(:)
    real(8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
    real(8), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
    real(8), allocatable :: dens(:,:,:), prs(:,:,:), psi(:,:,:)
    real(8) :: w_sum_ion, w_sum_elec

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    call set_random_seed(base_seed + 1299721*rank)
    if (rank == 0) print *, '----------------- Start 2D XZ initialization -----------------'
    if (rank == 0) print *, 'Grid = ', nx, ' x ', ny, ' x ', nz
    if (rank == 0) print *, 'nIon / nElec = ', nIon, nElec

    base_ion = nIon / nprocs
    extra_ion = mod(nIon, nprocs)
    nIon_local = base_ion
    if (rank < extra_ion) nIon_local = nIon_local + 1
    ion_offset = rank * base_ion + min(rank, extra_ion)

    base_elec = nElec / nprocs
    extra_elec = mod(nElec, nprocs)
    nElec_local = base_elec
    if (rank < extra_elec) nElec_local = nElec_local + 1
    elec_offset = rank * base_elec + min(rank, extra_elec)

    allocate(Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz))
    allocate(Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz))
    allocate(dens(nx,ny,nz), prs(nx,ny,nz), psi(nx,ny,nz))
    allocate(ion_local(nIon_local), elec_local(nElec_local))

    if (rank == 0) call load_equilibrium(Bx, By, Bz, dens, prs, Ex, Ey, Ez, psi, ierr)
    call MPI_Bcast(Bx, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(By, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Bz, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ex, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ey, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ez, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dens, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(prs, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(psi, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    call init_species(ion_local, nIon_local, nIon, M_ION_TEST(), 1000.0d0, q, &
                      Bx, By, Bz, Ex, Ey, Ez, dens, prs, 'ion', w_sum_ion, rank, ierr)

    call set_random_seed(base_seed + 99999 + 1299721*rank)
    call init_species(elec_local, nElec_local, nElec, M_ELEC_TEST(), 1000.0d0, -q, &
                      Bx, By, Bz, Ex, Ey, Ez, dens, prs, 'elec', w_sum_elec, rank, ierr)

    if (rank == 0) print *, '[charge neutrality] Q_total =', q*(w_sum_ion - w_sum_elec)
    if (rank == 0) print *, '----------------- Start writing particle data -----------------'
    call write_particle_data(ion_local, elec_local, nIon_local, nElec_local, rank, nprocs, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, '----------------- End 2D XZ initialization -------------------'

    deallocate(ion_local, elec_local)
    deallocate(Bx, By, Bz, Ex, Ey, Ez, dens, prs, psi)
    call MPI_Finalize(ierr)

contains
    real(8) function M_ELEC_TEST()
        M_ELEC_TEST = 9.1093837139d-31
    end function M_ELEC_TEST

    real(8) function M_ION_TEST()
        M_ION_TEST = 25.0d0 * M_ELEC_TEST()
    end function M_ION_TEST

    subroutine set_random_seed(seed_base)
        integer, intent(in) :: seed_base
        integer :: j, local_seed_size
        integer, allocatable :: local_seed(:)
        call random_seed(size=local_seed_size)
        allocate(local_seed(local_seed_size))
        do j = 1, local_seed_size
            local_seed(j) = int(mod(dble(seed_base) + 97.0d0*dble(j), 2147483647.0d0))
        end do
        call random_seed(put=local_seed)
        deallocate(local_seed)
    end subroutine set_random_seed

    subroutine load_equilibrium(Bx, By, Bz, dens, prs, Ex, Ey, Ez, psi, ierr)
        real(8), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        real(8), intent(out) :: dens(nx,ny,nz), prs(nx,ny,nz), psi(nx,ny,nz)
        real(8), intent(out) :: Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz)
        integer, intent(inout) :: ierr

        call read_xz('out_xz/Bx_xz.txt', Bx, ierr)
        call read_xz('out_xz/By_xz.txt', By, ierr)
        call read_xz('out_xz/Bz_xz.txt', Bz, ierr)
        call read_xz('out_xz/Ex_xz.txt', Ex, ierr)
        call read_xz('out_xz/Ey_xz.txt', Ey, ierr)
        call read_xz('out_xz/Ez_xz.txt', Ez, ierr)
        call read_xz('out_xz/density_xz.txt', dens, ierr)
        call read_xz('out_xz/pressure_xz.txt', prs, ierr)
        call read_xz('out_xz/psi_xz.txt', psi, ierr)
    end subroutine load_equilibrium

    subroutine read_xz(filename, field, ierr)
        character(len=*), intent(in) :: filename
        real(8), intent(out) :: field(nx,ny,nz)
        integer, intent(inout) :: ierr
        integer :: unit_no, ios
        open(newunit=unit_no, file=filename, form='formatted', status='old', &
             action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR: cannot open ', trim(filename)
            call MPI_Abort(MPI_COMM_WORLD, 71, ierr)
        end if
        read(unit_no, *, iostat=ios) field
        close(unit_no)
        if (ios /= 0) then
            print *, 'ERROR: cannot read ', trim(filename)
            call MPI_Abort(MPI_COMM_WORLD, 72, ierr)
        end if
    end subroutine read_xz

    subroutine init_species(part, n_local, n_global, mass, temp, charge, &
                             Bx, By, Bz, Ex, Ey, Ez, dens, prs, species, &
                             weight_sum_global, rank, ierr)
        type(particle), intent(out) :: part(:)
        integer, intent(in) :: n_local, n_global, rank
        real(8), intent(in) :: mass, temp, charge
        real(8), intent(in) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        real(8), intent(in) :: Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz)
        real(8), intent(in) :: dens(nx,ny,nz), prs(nx,ny,nz)
        character(len=*), intent(in) :: species
        real(8), intent(out) :: weight_sum_global
        integer, intent(inout) :: ierr

        integer :: i, ix, iz, ix_idx, iz_idx, ix_m, ix_p, iz_m, iz_p
        real(8) :: xx, zz, accept, x_phys, z_phys, alpha, gamma
        real(8) :: dens_max, dens_interpolated, n_real, weight_const
        real(8) :: vol_cell, v3(3), Bxi, Byi, Bzi, Exi, Eyi, Ezi
        real(8) :: B_sq, grad_Px, grad_Py, grad_Pz, pressure_fraction
        real(8) :: vE_x, vE_y, vE_z, vDia_x, vDia_y, vDia_z
        real(8) :: drift_mag, v_th, limit_factor
        real(8) :: weight_sum_local, weighting
        real(8), allocatable :: n_chk(:,:,:), n_chk_cell(:,:,:)
        real(8) :: err_threshold, dens_floor, err_max, err_sum_sq, rmse_val
        integer :: n_cells_checked

        n_real = 0.0d0
        if (rank == 0) then
            vol_cell = dx * dy * dz
            do iz = 1, nz - 1
                do ix = 1, nx - 1
                    n_real = n_real + dens(ix,1,iz) * vol_cell
                end do
            end do
            print *, '[Background field based] ', trim(species), ' = ', n_real
        end if
        call MPI_Bcast(n_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (n_global <= 0) call MPI_Abort(MPI_COMM_WORLD, 73, ierr)
        weight_const = n_real / dble(n_global)
        dens_max = maxval(dens)

        do i = 1, n_local
            do
                call random_number(xx); call random_number(zz); call random_number(accept)
                xx = min(xx, 1.0d0 - 1.0d-14)
                zz = min(zz, 1.0d0 - 1.0d-14)
                x_phys = x_min + xx * (x_max - x_min)
                z_phys = z_min + zz * (z_max - z_min)
                ix_idx = max(1, min(nx-1, int((x_phys-x_min)/dx) + 1))
                iz_idx = max(1, min(nz-1, int((z_phys-z_min)/dz) + 1))
                alpha = (x_phys - (x_min + dble(ix_idx-1)*dx)) / dx
                gamma = (z_phys - (z_min + dble(iz_idx-1)*dz)) / dz
                dens_interpolated = bilinear(dens, ix_idx, iz_idx, alpha, gamma)
                if (accept <= dens_interpolated / dens_max) exit
            end do

            part(i)%x = x_phys
            part(i)%y = y_plane
            part(i)%z = z_phys
            v3 = sample_maxwell(temp, mass)
            Bxi = bilinear(Bx, ix_idx, iz_idx, alpha, gamma)
            Byi = bilinear(By, ix_idx, iz_idx, alpha, gamma)
            Bzi = bilinear(Bz, ix_idx, iz_idx, alpha, gamma)
            Exi = bilinear(Ex, ix_idx, iz_idx, alpha, gamma)
            Eyi = bilinear(Ey, ix_idx, iz_idx, alpha, gamma)
            Ezi = bilinear(Ez, ix_idx, iz_idx, alpha, gamma)
            B_sq = Bxi*Bxi + Byi*Byi + Bzi*Bzi + 1.0d-20

            pressure_fraction = 0.5d0
            ix_m = max(1, ix_idx-1); ix_p = min(nx, ix_idx+1)
            iz_m = max(1, iz_idx-1); iz_p = min(nz, iz_idx+1)
            grad_Px = pressure_fraction * (prs(ix_p,1,iz_idx) - prs(ix_m,1,iz_idx)) / &
                      (dble(ix_p-ix_m) * dx)
            grad_Py = 0.0d0
            grad_Pz = pressure_fraction * (prs(ix_idx,1,iz_p) - prs(ix_idx,1,iz_m)) / &
                      (dble(iz_p-iz_m) * dz)

            vE_x = (Eyi*Bzi - Ezi*Byi) / B_sq
            vE_y = (Ezi*Bxi - Exi*Bzi) / B_sq
            vE_z = (Exi*Byi - Eyi*Bxi) / B_sq
            vDia_x = (Byi*grad_Pz - Bzi*grad_Py) / (charge*dens_interpolated*B_sq)
            vDia_y = (Bzi*grad_Px - Bxi*grad_Pz) / (charge*dens_interpolated*B_sq)
            vDia_z = (Bxi*grad_Py - Byi*grad_Px) / (charge*dens_interpolated*B_sq)

            drift_mag = sqrt((vE_x+vDia_x)**2 + (vE_y+vDia_y)**2 + (vE_z+vDia_z)**2)
            v_th = sqrt(2.0d0 * temp * q / mass)
            limit_factor = 1.0d0
            if (drift_mag > 2.0d0*v_th) limit_factor = 2.0d0*v_th/drift_mag
            part(i)%vx = v3(1) + (vE_x+vDia_x)*limit_factor
            part(i)%vy = v3(2) + (vE_y+vDia_y)*limit_factor
            part(i)%vz = v3(3) + (vE_z+vDia_z)*limit_factor
            part(i)%mass = mass
            part(i)%charge = charge
            part(i)%weight = weight_const
        end do

        weight_sum_local = 0.0d0
        do i = 1, n_local
            weight_sum_local = weight_sum_local + part(i)%weight
        end do
        call MPI_Reduce(weight_sum_local, weight_sum_global, 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        weighting = 0.0d0
        if (n_local > 0) weighting = part(1)%weight
        if (rank == 0) print *, '[Hyperparticle Weighting] ', trim(species), ' = ', weighting

        allocate(n_chk(nx,ny,nz), n_chk_cell(nx,ny,nz))
        n_chk = 0.0d0; n_chk_cell = 0.0d0
        call deposit_weights(part, n_local, n_chk)
        call MPI_Allreduce(MPI_IN_PLACE, n_chk, nx*ny*nz, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        if (rank == 0) then
            vol_cell = dx * dy * dz
            do iz = 1, nz-1
                do ix = 1, nx-1
                    n_chk_cell(ix,1,iz) = 0.25d0 * &
                        (n_chk(ix,1,iz) + n_chk(ix+1,1,iz) + &
                         n_chk(ix,1,iz+1) + n_chk(ix+1,1,iz+1)) / vol_cell
                end do
            end do
            dens_floor = 0.1d0 * dens_max
            err_max = 0.0d0; err_sum_sq = 0.0d0; n_cells_checked = 0
            do iz = 1, nz-1
                do ix = 1, nx-1
                    if (dens(ix,1,iz) > dens_floor) then
                        err_threshold = abs(n_chk_cell(ix,1,iz)-dens(ix,1,iz)) / dens(ix,1,iz)
                        err_max = max(err_max, err_threshold)
                        err_sum_sq = err_sum_sq + (n_chk_cell(ix,1,iz)-dens(ix,1,iz))**2
                        n_cells_checked = n_cells_checked + 1
                    end if
                end do
            end do
            rmse_val = 0.0d0
            if (n_cells_checked > 0) rmse_val = sqrt(err_sum_sq/dble(n_cells_checked))
            print *, '[', trim(species), '] init inaccuracies = ', err_max
            print *, '[', trim(species), '] RMSE of density = ', rmse_val
        end if
        deallocate(n_chk, n_chk_cell)
    end subroutine init_species

    real(8) function bilinear(field, ix, iz, alpha, gamma)
        real(8), intent(in) :: field(nx,ny,nz), alpha, gamma
        integer, intent(in) :: ix, iz
        bilinear = field(ix,1,iz)*(1.0d0-alpha)*(1.0d0-gamma) + &
                   field(ix+1,1,iz)*alpha*(1.0d0-gamma) + &
                   field(ix,1,iz+1)*(1.0d0-alpha)*gamma + &
                   field(ix+1,1,iz+1)*alpha*gamma
    end function bilinear

    subroutine deposit_weights(p, npart, w_grid)
        type(particle), intent(in) :: p(:)
        integer, intent(in) :: npart
        real(8), intent(inout) :: w_grid(nx,ny,nz)
        integer :: k, ix, iz
        real(8) :: alpha, gamma, x_i0, z_k0

        do k = 1, npart
            ix = max(1, min(nx-1, int((p(k)%x-x_min)/dx)+1))
            iz = max(1, min(nz-1, int((p(k)%z-z_min)/dz)+1))
            x_i0 = x_min + dble(ix-1)*dx
            z_k0 = z_min + dble(iz-1)*dz
            alpha = max(0.0d0, min(1.0d0, (p(k)%x-x_i0)/dx))
            gamma = max(0.0d0, min(1.0d0, (p(k)%z-z_k0)/dz))
            w_grid(ix,1,iz) = w_grid(ix,1,iz) + p(k)%weight*(1.0d0-alpha)*(1.0d0-gamma)
            w_grid(ix+1,1,iz) = w_grid(ix+1,1,iz) + p(k)%weight*alpha*(1.0d0-gamma)
            w_grid(ix,1,iz+1) = w_grid(ix,1,iz+1) + p(k)%weight*(1.0d0-alpha)*gamma
            w_grid(ix+1,1,iz+1) = w_grid(ix+1,1,iz+1) + p(k)%weight*alpha*gamma
        end do
    end subroutine deposit_weights

    subroutine write_particle_data(ion_local, elec_local, n_ion_local, n_elec_local, &
                                   rank, nprocs, ierr)
        type(particle), intent(in) :: ion_local(:), elec_local(:)
        integer, intent(in) :: n_ion_local, n_elec_local, rank, nprocs
        integer, intent(inout) :: ierr
        character(len=64) :: species
        character(len=10) :: open_status
        integer :: p, species_index

        do p = 0, nprocs-1
            if (rank == p) then
                open_status = 'append'
                if (p == 0) open_status = 'replace'
                do species_index = 1, 2
                    if (species_index == 1) then
                        species = 'ion'
                        call write_species_components(species, ion_local, n_ion_local, open_status, rank, ierr)
                    else
                        species = 'elec'
                        call write_species_components(species, elec_local, n_elec_local, open_status, rank, ierr)
                    end if
                end do
            end if
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do
    end subroutine write_particle_data

    subroutine write_species_components(species, p_array, n_local, open_status, rank, ierr)
        character(len=*), intent(in) :: species, open_status
        type(particle), intent(in) :: p_array(:)
        integer, intent(in) :: n_local, rank
        integer, intent(inout) :: ierr
        character(len=128) :: filename
        character(len=6), parameter :: components(6) = [character(len=6) :: &
            'x', 'z', 'vx', 'vy', 'vz', 'weight']
        integer :: c

        do c = 1, size(components)
            write(filename,'("init_particles/",A,"_",A,".bin")') trim(species), trim(components(c))
            call write_component_binary(filename, p_array, n_local, components(c), &
                                        open_status, rank, MPI_COMM_WORLD, ierr)
        end do
    end subroutine write_species_components

    subroutine write_component_binary(filename, p_array, n_local, component_name, &
                                      open_status, rank, comm, ierr)
        character(len=*), intent(in) :: filename, component_name, open_status
        type(particle), intent(in) :: p_array(:)
        integer, intent(in) :: n_local, rank, comm
        integer, intent(inout) :: ierr
        integer, parameter :: file_unit = 30
        integer :: i, io_stat, alloc_stat
        real(8), allocatable :: data_buffer(:)

        if (n_local == 0) return
        allocate(data_buffer(n_local), stat=alloc_stat)
        if (alloc_stat /= 0) then
            print *, 'ERROR: allocation failed for ', trim(filename), ' on rank ', rank
            call MPI_Abort(comm, 81, ierr)
        end if
        do i = 1, n_local
            select case (trim(component_name))
            case ('x'); data_buffer(i) = p_array(i)%x
            case ('z'); data_buffer(i) = p_array(i)%z
            case ('vx'); data_buffer(i) = p_array(i)%vx
            case ('vy'); data_buffer(i) = p_array(i)%vy
            case ('vz'); data_buffer(i) = p_array(i)%vz
            case ('weight'); data_buffer(i) = p_array(i)%weight
            end select
        end do
        if (trim(open_status) == 'append') then
            open(unit=file_unit, file=trim(filename), status='old', form='unformatted', &
                 access='stream', position='append', action='write', iostat=io_stat)
        else
            open(unit=file_unit, file=trim(filename), status='replace', form='unformatted', &
                 access='stream', action='write', iostat=io_stat)
        end if
        if (io_stat /= 0) then
            print *, 'ERROR: cannot open ', trim(filename), ' for writing.'
            call MPI_Abort(comm, 82, ierr)
        end if
        write(file_unit, iostat=io_stat) data_buffer
        close(file_unit)
        if (io_stat /= 0) then
            print *, 'ERROR: cannot write ', trim(filename)
            call MPI_Abort(comm, 83, ierr)
        end if
        deallocate(data_buffer)
    end subroutine write_component_binary
end program pic_equilibrium
