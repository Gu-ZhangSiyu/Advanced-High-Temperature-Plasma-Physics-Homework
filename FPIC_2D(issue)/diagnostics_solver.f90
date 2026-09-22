!***************************************************************
!   PIC诊断与场/粒子数据输出子程序
!
!	调用示例
!   call Output_Controller(idyn=idyn, dt=dt, freq_ion=10, freq_elec=100, freq_field=200, sample_rate_ion=10, sample_rate_elec=100, &
!                           out_ion=.true., out_elec=.true., out_field=.true., ion_local=ion_local, nIon_local=nIon_local, &
!                           elec_local=elec_local, nElec_local=nElec_local, Bx=Bx, By=By, Bz=Bz, Ex=Ex, Ey=Ey, Ez=Ez, &
!                           domain=domain)
!
!	参数设置：
!   freq_ion / freq_elec=10 / freq_field = 100	离子 / 电子 / 场 输出频率
!   sample_rate_ion / sample_rate_elec			粒子降采样率(每隔多少个ID记录一个粒子)
!   out_ion / out_elec / out_field =.true.		输出离子 / 电子 / 场 数据
!   
!   传递值说明:
!   idyn	当前的主循环时间步数	标量值，例如 1000
!   dt		物理时间步长(s)		标量值，例如 1.0d-12
!   freq_ion/freq_elec			离子/电子位置与速度的输出频率步数			标量值，例如 500 (每500步输出一次)
!   freq_field					电磁场三维矩阵的输出频率步数				标量值，例如 200 (每200步输出一次)
!   sample_rate_ion/sample_rate_elec	离子/电子空间采样率(按ID过滤)	标量值，100 表示仅记录ID能被100整除的粒子
!   out_ion/out_elec/out_field			是否开启离子/电子/场数据输出		布尔值，.true. 开启，.false. 关闭
!   ion_local/elec_local				本地进程的离子/电子群组			派生类型数组
!   nIon_local/nElec_local				本地进程的离子/电子总数			标量值
!   Bx,By,Bz,Ex,Ey,Ez					磁场/电场分量					包含卫星层的三维数组 (0:nx+1, 0:ny+1, 0:nz+1)
!***************************************************************
module diagnostics_solver
    use mpi_f08; use mesh_parameter; use physicals_parameter; use pic_types
    use Domain_Decomposition_Module, only: domain_decomposition
    implicit none

contains

	subroutine Output_Controller(idyn, dt, freq_ion, freq_elec, freq_field, &
                                sample_rate_ion, sample_rate_elec, &
                                out_ion, out_elec, out_field, &
                                ion_local, nIon_local, elec_local, nElec_local, &
                                Bx, By, Bz, Ex, Ey, Ez, domain)
        implicit none
        logical, intent(in) :: out_ion, out_elec, out_field
        type(particle), intent(in) :: ion_local(:), elec_local(:)
        integer, intent(in) :: idyn, freq_ion, freq_elec, freq_field, sample_rate_ion, sample_rate_elec, nIon_local, nElec_local
        real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:), Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:), dt
        type(domain_decomposition), intent(in) :: domain
        real(8) :: avg_r, std_r, ke, wtot_local

        if (out_ion .and. (mod(idyn, freq_ion) == 0 .or. idyn == 1)) then
            call compute_statistics(ion_local, nIon_local, avg_r, std_r, ke, wtot_local, idyn, dt, 'ion')
            call dump_dynamic_binary(ion_local, nIon_local, 'ions', idyn, sample_rate_ion)
        end if
        if (out_elec .and. (mod(idyn, freq_elec) == 0 .or. idyn == 1)) then
            call compute_statistics(elec_local, nElec_local, avg_r, std_r, ke, wtot_local, idyn, dt, 'elec')
            call dump_dynamic_binary(elec_local, nElec_local, 'electrons', idyn, sample_rate_elec)
        end if
        if (out_field .and. (mod(idyn, freq_field) == 0 .or. idyn == 1)) then
            call dump_electromagnetic_fields_text(Bx, By, Bz, Ex, Ey, Ez, idyn, domain)
        end if
    end subroutine Output_Controller
!*********************************************************************
subroutine compute_statistics(part, nP_local, avg_r, std_r, ke_tot, wtot_local, idyn, dt, species)	! 能量-粒子数检查和输出
    use mpi_f08; use physicals_parameter; use pic_types
    implicit none
	type(particle), intent(in) :: part(:)
	integer, intent(in) :: nP_local, idyn
    integer :: i, rank, ierr
    real(8), intent(in) :: dt
	real(8), intent(out):: avg_r, std_r, ke_tot, wtot_local
	real(8) :: w, delta, rad_pos, y_ke, t_ke, v2, ke_global, wtot_global, current_time, species_mass, &
				mean, m2, wtot, ke_sum, c_ke
    character(len=*), intent(in) :: species
	logical, save :: is_ion_initialized = .false.
    logical, save :: is_elec_initialized = .false.
	integer, save :: u_ke_ion, u_n_ion, u_ke_elec, u_n_elec

    mean = 0d0; m2 = 0d0; wtot = 0d0; ke_sum = 0d0; c_ke = 0d0
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (trim(species) == 'ion') then
        species_mass = M_ION
    else
        species_mass = M_ELEC
    end if

	do i = 1, nP_local
		w = part(i)%weight
		if (w <= 0d0) cycle
        rad_pos = sqrt(part(i)%x**2 + part(i)%z**2)		! XZ平面内的径向位置

		! Welford 在线更新均值/方差 (针对R)
		delta = rad_pos - mean
        wtot = wtot + w
        mean  = mean + w * delta / wtot
        m2 = m2  + w * delta * (rad_pos - mean)

		! 动能 Kahan 补偿求和 (3D速度)
		v2 = part(i)%vx**2 + part(i)%vy**2 + part(i)%vz**2
        y_ke = 0.5d0 * species_mass * v2 * w - c_ke
        t_ke = ke_sum + y_ke
        c_ke = (t_ke - ke_sum) - y_ke
        ke_sum = t_ke
	end do

	if (wtot <= 0d0) then
		avg_r = 0d0; std_r = 0d0; ke_tot = 0d0
	else
		avg_r = mean
		std_r = sqrt(max(m2 / wtot, 0d0))
		ke_tot = ke_sum
	end if
	wtot_local = wtot

	call MPI_Allreduce(ke_tot, ke_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
	call MPI_Allreduce(wtot_local, wtot_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        current_time = real(idyn, 8) * dt * 1.0d12	! 转化为皮秒
		if (trim(species) == 'ion') then
            if (.not. is_ion_initialized) then
                open(newunit=u_ke_ion, file='diagnostics/KE_ions.txt', status='unknown', action='write', position='append')
                open(newunit=u_n_ion, file='diagnostics/N_ions.txt',  status='unknown', action='write', position='append')
                is_ion_initialized = .true.
            end if
            write(u_ke_ion, '(E22.15, 1X, E22.15)') current_time, ke_global; flush(u_ke_ion)
            write(u_n_ion, '(E22.15, 1X, E22.15)') current_time, wtot_global; flush(u_n_ion)
        else if (trim(species) == 'elec') then
            if (.not. is_elec_initialized) then
                open(newunit=u_ke_elec, file='diagnostics/KE_electrons.txt', status='unknown', action='write', position='append')
                open(newunit=u_n_elec, file='diagnostics/N_electrons.txt',  status='unknown', action='write', position='append')
                is_elec_initialized = .true.
            end if
            write(u_ke_elec, '(E22.15, 1X, E22.15)') current_time, ke_global; flush(u_ke_elec)
            write(u_n_elec, '(E22.15, 1X, E22.15)') current_time, wtot_global; flush(u_n_elec)
        end if
    end if
end subroutine compute_statistics
!*********************************************************************
subroutine dump_dynamic_binary(part, nLocal, species_in, step, sample_rate)	! 粒子位置与速度输出
    use mpi_f08; use mesh_parameter; use pic_types
    implicit none
    type(particle), intent(in) :: part(:)
	integer :: ierr, rank, i, active_p_count
	type(MPI_File) :: fh_pos, fh_vel
	type(MPI_Status) :: status
	integer, intent(in) :: nLocal, step, sample_rate
    integer(kind=8) :: local_count, global_offset
    integer(kind=MPI_OFFSET_KIND) :: disp
	integer, allocatable :: valid_idx(:)
    character(len=64) :: fname_pos, fname_vel
    character(len=12) :: species
	character(len=*), intent(in):: species_in
    real(kind=8), allocatable :: data_pos(:,:), data_vel(:,:)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (sample_rate < 1) call MPI_Abort(MPI_COMM_WORLD, 42, ierr)
    species = adjustl(species_in)
    write(fname_pos, '("dynamic_",A,"_position/",I0,".bin")') trim(species), step
    write(fname_vel, '("dynamic_",A,"_velocity/",I0,".bin")') trim(species), step

	! 根据权重和ID采样建立输出索引；权重本身不再写入动态输出。
	allocate(valid_idx(nLocal))
    active_p_count = 0
    do i = 1, nLocal
        if (part(i)%weight /= 0.0d0 .and. mod(part(i)%id, int(sample_rate, 8)) == 0_8) then
            active_p_count = active_p_count + 1; valid_idx(active_p_count) = i  
        end if
    end do
    local_count = int(active_p_count, kind=8)

    ! 确定全局偏移量
    call MPI_Exscan(local_count, global_offset, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rank == 0) global_offset = 0_8
    disp = global_offset * 32_8

    ! 数据打包
	if (active_p_count > 0) then
        allocate(data_pos(4, active_p_count)); allocate(data_vel(4, active_p_count))
        do i = 1, active_p_count
            data_pos(1, i) = dble(part(valid_idx(i))%id)
            data_pos(2, i) = part(valid_idx(i))%x
            data_pos(3, i) = part(valid_idx(i))%y
            data_pos(4, i) = part(valid_idx(i))%z
            data_vel(1, i) = dble(part(valid_idx(i))%id)
            data_vel(2, i) = part(valid_idx(i))%vx
            data_vel(3, i) = part(valid_idx(i))%vy
            data_vel(4, i) = part(valid_idx(i))%vz
        end do
    else
        ! MPI-IO 仍要求提供有效的缓冲地址；count 保持为 0，dummy 不会写入文件。
        allocate(data_pos(4, 1), data_vel(4, 1))
        data_pos = 0.0d0; data_vel = 0.0d0
    end if
    deallocate(valid_idx)

    ! MPI-IO集体写入
    call MPI_File_open(MPI_COMM_WORLD, trim(fname_pos), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh_pos, ierr)
    call MPI_File_set_size(fh_pos, 0_MPI_OFFSET_KIND, ierr)
    call MPI_File_write_at_all(fh_pos, disp, data_pos, active_p_count * 4, MPI_DOUBLE_PRECISION, status, ierr)
    call MPI_File_close(fh_pos, ierr)
    call MPI_File_open(MPI_COMM_WORLD, trim(fname_vel), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh_vel, ierr)
    call MPI_File_set_size(fh_vel, 0_MPI_OFFSET_KIND, ierr)
    call MPI_File_write_at_all(fh_vel, disp, data_vel, active_p_count * 4, MPI_DOUBLE_PRECISION, status, ierr)
    call MPI_File_close(fh_vel, ierr)
    deallocate(data_pos, data_vel)
end subroutine dump_dynamic_binary
!*********************************************************************
subroutine dump_electromagnetic_fields_text(Bx, By, Bz, Ex, Ey, Ez, step, domain) ! 129x129 XZ场文本输出
    implicit none
    integer, intent(in) :: step
    real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:)
    real(8), intent(in) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    character(len=128) :: filename

    write(filename, '("dynamic_fields_Bx/",I0,".txt")') step
    call write_field_matrix_text(filename, Bx, domain)
    write(filename, '("dynamic_fields_By/",I0,".txt")') step
    call write_field_matrix_text(filename, By, domain)
    write(filename, '("dynamic_fields_Bz/",I0,".txt")') step
    call write_field_matrix_text(filename, Bz, domain)
    write(filename, '("dynamic_fields_Ex/",I0,".txt")') step
    call write_field_matrix_text(filename, Ex, domain)
    write(filename, '("dynamic_fields_Ey/",I0,".txt")') step
    call write_field_matrix_text(filename, Ey, domain)
    write(filename, '("dynamic_fields_Ez/",I0,".txt")') step
    call write_field_matrix_text(filename, Ez, domain)
end subroutine dump_electromagnetic_fields_text
!*********************************************************************
subroutine write_field_matrix_text(filename, field, domain) ! 每个文件严格为nz行、nx列
    implicit none
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: field(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    integer :: rank, nprocs, writer, i, k, unit_no, ios, ierr

    call MPI_Comm_rank(domain%comm, rank, ierr)
    call MPI_Comm_size(domain%comm, nprocs, ierr)
    if (rank == 0) then
        open(newunit=unit_no, file=trim(filename), form='formatted', status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) call MPI_Abort(domain%comm, 45, ierr)
        close(unit_no)
    end if
    call MPI_Barrier(domain%comm, ierr)

    ! 行顺序为全局 z=1...nz，列顺序为 x=1...nx；不写表头，便于直接读成129x129矩阵。
    do writer = 0, nprocs - 1
        if (rank == writer) then
            open(newunit=unit_no, file=trim(filename), form='formatted', status='old', &
                 position='append', action='write', iostat=ios)
            if (ios /= 0) call MPI_Abort(domain%comm, 45, ierr)
            do k = 1, domain%nz_local
                write(unit_no, '(*(ES24.16,1X))', iostat=ios) (field(i,1,k), i=1,nx)
                if (ios /= 0) call MPI_Abort(domain%comm, 45, ierr)
            end do
            close(unit_no)
        end if
        call MPI_Barrier(domain%comm, ierr)
    end do
end subroutine write_field_matrix_text
end module diagnostics_solver
