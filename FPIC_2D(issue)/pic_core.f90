!***************************************************************
!   Full-PIC 主程序核心计算
!*********************************************************************
program pic_equilibrium
    use mpi_f08
	use mesh_parameter; use physicals_parameter; use pic_types; use control_parameter
    use particle_solver_yee, only: Push_Particle
    use diagnostics_solver, only: Output_Controller
	use data_loader_module, only: load_equilibrium, read_particle_data
	use restart_manager_module_CARA; use beam_injection_module
	use Particle_Operations; use Utility_Functions; use Ghost_Cells_Module
    use Domain_Decomposition_Module
    use Particle_Migration_Module
    use Load_Balance_Module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
	integer :: rank, nprocs, ierr, idyn, max_ion_buffer_local, max_elec_buffer_local, &
			nIon_local, nElec_local, base_Ion, extra_Ion, base_Elec, extra_Elec, ion_offset, elec_offset
	type(domain_decomposition) :: domain
	integer(kind=8) :: global_nIon_alive, global_nElec_alive, global_inj_ion, global_inj_elec
	real(8), allocatable :: Bx0(:,:,:), By0(:,:,:), Bz0(:,:,:), Ex0(:,:,:), Ey0(:,:,:), Ez0(:,:,:)
    real(8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:), Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), Jx(:,:,:), Jy(:,:,:), Jz(:,:,:), &
							Bx_half(:,:,:), By_half(:,:,:), Bz_half(:,:,:), Bx_old(:,:,:), By_old(:,:,:), Bz_old(:,:,:)
	real(8) :: qm_ion, qm_elec, coef_curl, coef_current, cara_rf_time, cara_rf_Ey
    real(8) :: frac_ion = 0.0d0, frac_elec = 0.0d0
    type(particle), pointer :: ion_local(:), elec_local(:)
    real(8), allocatable :: Jx_ion(:,:,:), Jy_ion(:,:,:), Jz_ion(:,:,:)
    real(8) :: dt_ion, ion_inject_dt
    integer :: ion_block_end, ion_block_steps, last_balance_step, injection_step, ion_inject_steps
    logical :: ion_time_ready

	! 初始化MPI
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call Validate_Controls()
    call Initialize_Domain(domain, MPI_COMM_WORLD)
    if (enable_beam_injection) then
        if (injection_duration <= 0.0d0 .or. min(total_target_ions, total_target_elec) < 0.0d0 .or. &
            min(weight_ion, weight_elec) <= 0.0d0 .or. &
            Owner_Of_Position(z_min + inject_z_offset * dz, domain) == MPI_PROC_NULL) then
            if (rank == 0) print *, 'ERROR: invalid injection duration, count, weight or injection plane.'
            call MPI_Abort(MPI_COMM_WORLD, 51, ierr)
        end if
    end if
	if (rank == 0) print *, '----------------- Start Code -----------------'
	call Report_Domain_Decomposition(domain)
	
	! 每个rank仅分配本地z-slab加两个ghost平面
    allocate(Bx0(0:nx+1,0:ny+1,0:domain%nz_local+1), By0(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Bz0(0:nx+1,0:ny+1,0:domain%nz_local+1), Ex0(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Ey0(0:nx+1,0:ny+1,0:domain%nz_local+1), Ez0(0:nx+1,0:ny+1,0:domain%nz_local+1))
    allocate(Bx(0:nx+1,0:ny+1,0:domain%nz_local+1), By(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Bz(0:nx+1,0:ny+1,0:domain%nz_local+1), Ex(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Ey(0:nx+1,0:ny+1,0:domain%nz_local+1), Ez(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Bx_half(0:nx+1,0:ny+1,0:domain%nz_local+1), By_half(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Bz_half(0:nx+1,0:ny+1,0:domain%nz_local+1), Bx_old(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            By_old(0:nx+1,0:ny+1,0:domain%nz_local+1), Bz_old(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Jx(0:nx+1,0:ny+1,0:domain%nz_local+1), Jy(0:nx+1,0:ny+1,0:domain%nz_local+1), &
            Jz(0:nx+1,0:ny+1,0:domain%nz_local+1))
    Bx0 = 0.0d0; By0 = 0.0d0; Bz0 = 0.0d0; Ex0 = 0.0d0; Ey0 = 0.0d0; Ez0 = 0.0d0
    Bx = 0.0d0; By = 0.0d0; Bz = 0.0d0; Ex = 0.0d0; Ey = 0.0d0; Ez = 0.0d0; Jx = 0.0d0; Jy = 0.0d0; Jz = 0.0d0
    Bx_half = 0.0d0; By_half = 0.0d0; Bz_half = 0.0d0; Bx_old = 0.0d0; By_old = 0.0d0; Bz_old = 0.0d0

	! 重启的优先级高于后台加载, 先记录分区文件, 读取后按空间所有权重新分配（MPI大小发生变化时也是如此）
    global_nIon_alive = 0_8; global_nElec_alive = 0_8
    if (restart_step < 0 .or. nIon < 0_8 .or. nElec < 0_8) call MPI_Abort(MPI_COMM_WORLD, 40, ierr)
    if (restart_step > 0) then
        call get_active_particle_count('ions', restart_step, global_nIon_alive)
        call get_active_particle_count('electrons', restart_step, global_nElec_alive)
    else if (load_background_particles) then
        global_nIon_alive = nIon; global_nElec_alive = nElec
    end if

	! 需显式拒绝溢出
    if (max(global_nIon_alive, global_nElec_alive) > int(huge(0)/8,8)) then
        if (rank == 0) print *, 'ERROR: particle count exceeds current MPI packing count limit.'
        call MPI_Abort(MPI_COMM_WORLD, 41, ierr)
    end if
    base_Ion = int(global_nIon_alive) / nprocs; extra_Ion = mod(int(global_nIon_alive), nprocs); nIon_local = base_Ion
    if (rank < extra_Ion) nIon_local = nIon_local + 1
    ion_offset = rank * base_Ion + min(rank, extra_Ion); base_Elec = int(global_nElec_alive) / nprocs
    extra_Elec = mod(int(global_nElec_alive), nprocs); nElec_local = base_Elec
    if (rank < extra_Elec) nElec_local = nElec_local + 1
    elec_offset = rank * base_Elec + min(rank, extra_Elec)

    ! 记录分配空间动态扩展
    max_ion_buffer_local = max(1024, nIon_local); max_elec_buffer_local = max(1024, nElec_local)
    allocate(ion_local(max_ion_buffer_local), elec_local(max_elec_buffer_local))
    global_inj_ion = global_nIon_alive; global_inj_elec = global_nElec_alive
    frac_ion = 0.0d0; frac_elec = 0.0d0

	! 每个rank从全局文件读取自己的连续z-slab
    call load_equilibrium(Bx0, By0, Bz0, Ex0, Ey0, Ez0, domain)
    call Ghost_Cells(Bx0, 'Bx', domain); call Ghost_Cells(By0, 'By', domain); call Ghost_Cells(Bz0, 'Bz', domain)
    call Ghost_Cells(Ex0, 'Ex', domain); call Ghost_Cells(Ey0, 'Ey', domain); call Ghost_Cells(Ez0, 'Ez', domain)
	
    ! 初始化扰动场
	if (restart_step > 0) then		! 重启读取变化场
		call load_restart_fields(Bx, By, Bz, Ex, Ey, Ez, restart_step, domain)
    end if

    ! 初始化卫星层
    call Ghost_Cells(Bx, 'Bx', domain); call Ghost_Cells(By, 'By', domain); call Ghost_Cells(Bz, 'Bz', domain)
    call Ghost_Cells(Ex, 'Ex', domain); call Ghost_Cells(Ey, 'Ey', domain); call Ghost_Cells(Ez, 'Ez', domain)
	
    if (restart_step > 0) then
        call read_restart_particles(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, &
                                    rank, nprocs, ierr, restart_step, legacy_restart_weight_ion, legacy_restart_weight_elec)
        call initialize_injection_state(ion_local, nIon_local, global_inj_ion, frac_ion)
        call initialize_injection_state(elec_local, nElec_local, global_inj_elec, frac_elec)
    else if (load_background_particles) then
        call read_particle_data(ion_local, elec_local, nIon_local, ion_offset, nElec_local, elec_offset, rank, nprocs, ierr)
    end if
	
    ! 移除最初处于非活动状态或位于范围外的粒子
    call Remove_Outside_Particles(ion_local, nIon_local)
    call Remove_Outside_Particles(elec_local, nElec_local)
    call Redistribute_Particles_Global(ion_local, nIon_local, max_ion_buffer_local, domain)
    call Redistribute_Particles_Global(elec_local, nElec_local, max_elec_buffer_local, domain)
    call Report_Particle_Load(nIon_local, nElec_local, restart_step, domain)
	
	! 分数注入状态专属于注入平面所有者, 源值循环回小于1, 其他层级不得保留过期分数
    if (rank /= Owner_Of_Position(z_min + inject_z_offset * dz, domain)) then
		frac_ion = 0.0d0; frac_elec = 0.0d0
	end if
	if (sub_ratio > 1) then
        call Reset_Work_Field(Jx_ion); call Reset_Work_Field(Jy_ion); call Reset_Work_Field(Jz_ion)
    end if
    last_balance_step=restart_step
    if (enable_load_balance) call Rebalance(restart_step)
    ion_block_end=restart_step
    call Check_Subcycle_Restart()
	
    qm_ion = q / M_ION; qm_elec = -q / M_ELEC   ! 预计算荷质比
	coef_curl = dt * c2; coef_current = dt * inv_eps0

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (rank == 0) print *, '----------------- Start Main -----------------'
	if (rank == 0) then
        print *, 'Ion sub_ratio / dynamic load balance: ', sub_ratio, enable_load_balance
        print *, 'Load background / inject beam / RF: ', load_background_particles, enable_beam_injection, enable_RF
		print *, 'CARA RF frequency [MHz]: ', CARA_RF_omega / (2.0d0 * pi * 1.0d6)
		print *, 'CARA RF Ey peak [V/m]:  ', CARA_RF_Ey_amplitude
	end if

    ! 磁场时间与读取初始值同步
    if (restart_step == 0) call Update_B_FDTD(Bx, By, Bz, Ex, Ey, Ez, -dt/2.0d0, domain)

    ! 蛙跳主循环
    do idyn = restart_step + 1, nStepsDyn
        if (rank == 0 .and. mod(idyn, 100) == 0) print *, 'Step: ', idyn

        ! 当ie都达到相同的位置时间时才进行重划分
        if (enable_load_balance .and. idyn == ion_block_end+1) then
            if (idyn-1-last_balance_step >= load_balance_interval) then
                call Rebalance(idyn-1)
                last_balance_step=idyn-1
            end if
        end if

        ! 备份旧磁场 B(n-1/2)
        Bx_old = Bx; By_old = By; Bz_old = Bz

        ! 更新磁场 B(n-1/2) -> B(n+1/2)
        call Update_B_FDTD(Bx, By, Bz, Ex, Ey, Ez, dt, domain)

		! 构造时间中心总磁场
        Bx_half = 0.5d0 * (Bx + Bx_old) + Bx0; By_half = 0.5d0 * (By + By_old) + By0; Bz_half = 0.5d0 * (Bz + Bz_old) + Bz0
        Ex(:,:,:) = Ex + Ex0; Ey(:,:,:) = Ey + Ey0; Ez(:,:,:) = Ez + Ez0

		! 外加RF场, 计算均匀Ey(t), 每个rank对本地slab及ghost层叠加相同RF值
		cara_rf_Ey = 0.0d0
        if (enable_RF) then
			cara_rf_time = (dble(idyn) - 0.5d0) * dt
			cara_rf_Ey = CARA_RF_Ey_amplitude * sin(CARA_RF_omega * cara_rf_time + CARA_RF_phase)
			Ey(:,:,:) = Ey + cara_rf_Ey
        end if

        if (sub_ratio == 1) then
			! 束流连续同步注入
			if (enable_beam_injection) then
				call Inject_Beam_Species(ion_local, nIon_local, max_ion_buffer_local, dt, dt, idyn, domain, &
										M_ION, q, total_target_ions, weight_ion, E_beam_ion_keV, global_inj_ion, frac_ion)
				call Inject_Beam_Species(elec_local, nElec_local, max_elec_buffer_local, dt, dt, idyn, domain, &
										M_ELEC, -q, total_target_elec, weight_elec, E_beam_elec_keV, global_inj_elec, frac_elec)
			end if
		
			! 离子与电子推进
			call Run_Particle_Pusher(ion_local, nIon_local, Ex, Ey, Ez, Bx_half, By_half, Bz_half, dt, 0, qm_ion, domain)
			call Run_Particle_Pusher(elec_local, nElec_local, Ex, Ey, Ez, Bx_half, By_half, Bz_half, dt, 0, qm_elec, domain)

			! 跨界粒子迁移到新owner, 由新owner沉积
			call Migrate_Particles_Neighbors(ion_local, nIon_local, max_ion_buffer_local, domain, 'ion')
			call Migrate_Particles_Neighbors(elec_local, nElec_local, max_elec_buffer_local, domain, 'electron')
			
			! 电流沉积
			Jx = 0.0d0; Jy = 0.0d0; Jz = 0.0d0
			
			! 离子与电子电流叠加
			if (nIon_local > 0) call Deposit_Species_Yee(ion_local, nIon_local, Jx, Jy, Jz, q, dt, domain)
			if (nElec_local > 0) call Deposit_Species_Yee(elec_local, nElec_local, Jx, Jy, Jz, -q, dt, domain)
			call Remove_Outside_Particles(ion_local, nIon_local)
			call Remove_Outside_Particles(elec_local, nElec_local)
			call Reduce_And_Normalize_Current(Jx, Jy, Jz, dx, dy, dz, domain)
            ion_block_end=idyn
        else
            if (idyn == ion_block_end+1) then
                ion_block_steps=min(sub_ratio,nStepsDyn-idyn+1); ion_block_end=idyn+ion_block_steps-1
                dt_ion=dt*dble(ion_block_steps)
                if (enable_beam_injection) then
                    ion_inject_steps=0
                    do injection_step=idyn,ion_block_end
                        if (dble(injection_step)*dt > injection_duration) exit
                        ion_inject_steps=ion_inject_steps+1
                    end do
                    ion_inject_dt=dt*dble(ion_inject_steps)
                    if (ion_inject_dt > 0.0d0) then
                        call Inject_Beam_Species(ion_local,nIon_local,max_ion_buffer_local,ion_inject_dt,dt, &
                            idyn,domain,M_ION,q,total_target_ions,weight_ion,E_beam_ion_keV,global_inj_ion,frac_ion)
                    end if
                end if
                call Run_Particle_Pusher(ion_local,nIon_local,Ex,Ey,Ez,Bx_half,By_half,Bz_half, dt_ion,0,qm_ion,domain)
                Jx_ion=0.0d0; Jy_ion=0.0d0; Jz_ion=0.0d0
                call Deposit_Long_Orbits(ion_local,nIon_local,Jx_ion,Jy_ion,Jz_ion,q,dt_ion,domain)
                call Reduce_And_Normalize_Current(Jx_ion,Jy_ion,Jz_ion,dx,dy,dz,domain)
                call Remove_Outside_Particles(ion_local,nIon_local)
                call Redistribute_Particles_Global(ion_local,nIon_local,max_ion_buffer_local,domain)
            end if
            if (enable_beam_injection) then
                call Inject_Beam_Species(elec_local,nElec_local,max_elec_buffer_local,dt,dt,idyn,domain, &
                    M_ELEC,-q,total_target_elec,weight_elec,E_beam_elec_keV,global_inj_elec,frac_elec)
            end if
            call Run_Particle_Pusher(elec_local,nElec_local,Ex,Ey,Ez,Bx_half,By_half,Bz_half,dt,0,qm_elec,domain)
            call Migrate_Particles_Neighbors(elec_local,nElec_local,max_elec_buffer_local,domain,'electron')
            Jx=0.0d0; Jy=0.0d0; Jz=0.0d0
            if (nElec_local > 0) call Deposit_Species_Yee(elec_local,nElec_local,Jx,Jy,Jz,-q,dt,domain)
            call Remove_Outside_Particles(elec_local,nElec_local)
            call Reduce_And_Normalize_Current(Jx,Jy,Jz,dx,dy,dz,domain)
            Jx=Jx+Jx_ion; Jy=Jy+Jy_ion; Jz=Jz+Jz_ion
        end if
        ion_time_ready=(idyn == ion_block_end)

		! 扣除解析RF场和背景电场
		if (enable_RF) Ey(:,:,:) = Ey - cara_rf_Ey
        Ex(:,:,:) = Ex - Ex0; Ey(:,:,:) = Ey - Ey0; Ez(:,:,:) = Ez - Ez0
		
		! 总电流平滑
		call Apply_Binomial_Filter(Jx, Jy, Jz, domain)

        ! 更新电场 E(n) -> E(n+1)
        call Update_E_FDTD(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, domain)

        ! 诊断与输出				
		call Output_Controller(idyn=idyn, dt=dt, freq_ion=freq_ion, freq_elec=freq_elec, freq_field=freq_field, &
							sample_rate_ion=sample_rate_ion, sample_rate_elec=sample_rate_elec, &
							out_ion=out_ion .and. ion_time_ready, out_elec=out_elec, out_field=out_field, &
							ion_local=ion_local, nIon_local=nIon_local, &
							elec_local=elec_local, nElec_local=nElec_local, Bx=Bx, By=By, Bz=Bz, Ex=Ex, Ey=Ey, Ez=Ez, &
							domain=domain)
        call Write_Subcycle_State(idyn,ion_time_ready)
        if ((mod(idyn,1000) == 0 .or. idyn == 1) .and. ion_time_ready) call Report_Particle_Load(nIon_local,nElec_local,idyn,domain)
    end do
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (rank == 0) print *, '----------------- End Main -----------------'
	call MPI_Finalize(ierr)
contains
!*********************************************************************
! 分配/I/O前进行验证
subroutine Validate_Controls()
    if (sub_ratio < 1 .or. dt <= 0.0d0 .or. .not. ieee_is_finite(dt)) then
        if (rank == 0) print *, 'ERROR: sub_ratio must be >= 1 and dt positive/finite.'
        call MPI_Abort(MPI_COMM_WORLD,62,ierr)
    end if
    if (min(freq_ion,freq_elec,freq_field,sample_rate_ion,sample_rate_elec) < 1) then
        if (rank == 0) print *, 'ERROR: output frequencies and sampling rates must be positive.'
        call MPI_Abort(MPI_COMM_WORLD,62,ierr)
    end if
    if (sub_ratio > 1 .and. out_ion) then
        if (mod(freq_ion,sub_ratio) /= 0) then
            if (rank == 0) print *, 'ERROR: freq_ion must be a multiple of sub_ratio (ions output at macro boundaries).'
            call MPI_Abort(MPI_COMM_WORLD,62,ierr)
        end if
    end if
    if (enable_load_balance) then
        if (load_balance_interval < 1 .or. load_balance_threshold < 1.0d0 .or. &
            load_balance_min_gain < 0.0d0 .or. load_balance_min_gain >= 1.0d0 .or. &
            load_balance_cell_cost <= 0.0d0 .or. load_balance_electron_cost < 0.0d0 .or. &
            load_balance_ion_cost < 0.0d0 .or. load_balance_max_slab_factor < 1.0d0 .or. &
            .not. all(ieee_is_finite([load_balance_threshold,load_balance_min_gain,load_balance_cell_cost, &
                load_balance_electron_cost,load_balance_ion_cost,load_balance_max_slab_factor]))) then
            if (rank == 0) print *, 'ERROR: invalid load-balance controls.'
            call MPI_Abort(MPI_COMM_WORLD,62,ierr)
        end if
    end if
end subroutine Validate_Controls
!*********************************************************************
subroutine Reset_Work_Field(F)
    real(8), allocatable, intent(inout) :: F(:,:,:)
    if (allocated(F)) deallocate(F)
    allocate(F(0:nx+1,0:ny+1,0:domain%nz_local+1))
    F=0.0d0
end subroutine Reset_Work_Field
!*********************************************************************
subroutine Rebalance(step)
    integer, intent(in) :: step
    type(domain_decomposition) :: proposed, previous
    logical :: changed
    integer :: old_source,new_source,seed_size
    integer, allocatable :: seed(:)
    integer(kind=8) :: ids(2), all_ids(2)
    real(8) :: fractions(2), all_fractions(2), started
    started=MPI_Wtime()
    call Propose_Balanced_Domain(ion_local,nIon_local,elec_local,nElec_local,domain,proposed,changed,step)
    if (.not. changed) return
    previous=domain
    if (enable_beam_injection) then
        old_source=Owner_Of_Position(z_min+inject_z_offset*dz,previous)
        new_source=Owner_Of_Position(z_min+inject_z_offset*dz,proposed)
        ids=[global_inj_ion,global_inj_elec]; fractions=[frac_ion,frac_elec]
        call MPI_Allreduce(ids,all_ids,2,MPI_INTEGER8,MPI_MAX,domain%comm,ierr)
        call MPI_Allreduce(fractions,all_fractions,2,MPI_DOUBLE_PRECISION,MPI_SUM,domain%comm,ierr)
        global_inj_ion=all_ids(1); global_inj_elec=all_ids(2)
        frac_ion=0.0d0; frac_elec=0.0d0
        if (rank == new_source) then
            frac_ion=all_fractions(1); frac_elec=all_fractions(2)
        end if
        if (new_source /= old_source) then
            call random_seed(size=seed_size)
            allocate(seed(seed_size))
            if (rank == old_source) call random_seed(get=seed)
            call MPI_Bcast(seed,seed_size,MPI_INTEGER,old_source,domain%comm,ierr)
            if (rank == new_source) call random_seed(put=seed)
            deallocate(seed)
        end if
    end if
    call Remap_Field(Bx0,previous,proposed)
    call Remap_Field(By0,previous,proposed)
    call Remap_Field(Bz0,previous,proposed)
    call Remap_Field(Ex0,previous,proposed)
    call Remap_Field(Ey0,previous,proposed)
    call Remap_Field(Ez0,previous,proposed)
    call Remap_Field(Bx,previous,proposed)
    call Remap_Field(By,previous,proposed)
    call Remap_Field(Bz,previous,proposed)
    call Remap_Field(Ex,previous,proposed)
    call Remap_Field(Ey,previous,proposed)
    call Remap_Field(Ez,previous,proposed)
    domain=proposed
    call Ghost_Cells(Bx0,'Bx',domain)
    call Ghost_Cells(By0,'By',domain)
    call Ghost_Cells(Bz0,'Bz',domain)
    call Ghost_Cells(Ex0,'Ex',domain)
    call Ghost_Cells(Ey0,'Ey',domain)
    call Ghost_Cells(Ez0,'Ez',domain)
    call Ghost_Cells(Bx,'Bx',domain)
    call Ghost_Cells(By,'By',domain)
    call Ghost_Cells(Bz,'Bz',domain)
    call Ghost_Cells(Ex,'Ex',domain)
    call Ghost_Cells(Ey,'Ey',domain)
    call Ghost_Cells(Ez,'Ez',domain)
    call Reset_Work_Field(Bx_half)
    call Reset_Work_Field(By_half)
    call Reset_Work_Field(Bz_half)
    call Reset_Work_Field(Bx_old)
    call Reset_Work_Field(By_old)
    call Reset_Work_Field(Bz_old)
    call Reset_Work_Field(Jx)
    call Reset_Work_Field(Jy)
    call Reset_Work_Field(Jz)
    if (sub_ratio > 1) then
        call Reset_Work_Field(Jx_ion); call Reset_Work_Field(Jy_ion); call Reset_Work_Field(Jz_ion)
    end if
    call Redistribute_Particles_Global(ion_local,nIon_local,max_ion_buffer_local,domain)
    call Redistribute_Particles_Global(elec_local,nElec_local,max_elec_buffer_local,domain)
    call Report_Domain_Decomposition(domain)
    call Report_Particle_Load(nIon_local,nElec_local,step,domain)
    if (rank == 0) print *, 'balance transfer seconds=', MPI_Wtime()-started
end subroutine Rebalance
!*********************************************************************
! 辅助进程
subroutine Check_Subcycle_Restart()
    character(len=128) :: filename
    integer :: unit, ios, version, saved_step, saved_ratio, valid
    real(8) :: saved_dt
    logical :: exists
    if (restart_step == 0) return
    valid = 1
    if (rank == 0) then
        write(filename, '("for_restart/subcycle_",I0,".txt")') restart_step
        inquire(file=filename, exist=exists)
        if (exists) then
            open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            if (ios == 0) then
                read(unit, *, iostat=ios) version, saved_step, saved_ratio, saved_dt
                close(unit)
            end if
            if (ios /= 0) then
                valid = 0
            else if (version /= 1 .or. saved_step /= restart_step .or. saved_ratio /= sub_ratio .or. &
                     .not. ieee_is_finite(saved_dt) .or. saved_dt /= dt) then
                valid = 0
            end if
        else if (sub_ratio > 1) then
            valid = 0
        end if
        if (mod(restart_step, sub_ratio) /= 0) valid = 0
        if (valid == 0) print *, 'ERROR: restart requires a full macro boundary and matching subcycle_STEP.txt, ratio and dt.'
    end if
    call MPI_Bcast(valid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (valid == 0) call MPI_Abort(MPI_COMM_WORLD, 63, ierr)
end subroutine Check_Subcycle_Restart
!*********************************************************************
subroutine Write_Subcycle_State(step,ready)
    integer, intent(in) :: step
    logical, intent(in) :: ready
    character(len=128) :: filename
    integer :: unit,ios
    if (rank /= 0 .or. .not. out_field .or. .not. ready) return
    if (mod(step,sub_ratio) /= 0) return
    if (mod(step,freq_field) /= 0 .and. step /= 1) return
    write(filename,'("dynamic_fields_Ex/subcycle_",I0,".txt")') step
    open(newunit=unit,file=filename,status='replace',action='write',iostat=ios)
    if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD,64,ierr)
    write(unit,*,iostat=ios) 1,step,sub_ratio,dt
    if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD,64,ierr)
    close(unit)
end subroutine Write_Subcycle_State
!*********************************************************************
subroutine Update_B_FDTD(Bx, By, Bz, Ex, Ey, Ez, dt, domain)	! 磁场更新器
	use Ghost_Cells_Module
	use mesh_parameter, only: nx, nz, inv_dx, inv_dz
	use Domain_Decomposition_Module, only: domain_decomposition
	implicit none
	real(8), intent(inout) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:)
    real(8), intent(in) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
	real(8), intent(in) :: dt
	type(domain_decomposition), intent(in) :: domain
	real(8) :: dt_invdx, dt_invdz
	integer :: i, k, k_forward_end

	dt_invdx = dt * inv_dx; dt_invdz = dt * inv_dz
	k_forward_end = domain%nz_local
	if (domain%k_end == nz) k_forward_end = domain%nz_local - 1
	! Faraday 方程在 XZ 平面上的形式：
	! dBx/dt= dEy/dz, dBy/dt=-dEx/dz+dEz/dx, dBz/dt=-dEy/dx。
	do k = 1, k_forward_end
		do i = 1, nx
			Bx(i,1,k) = Bx(i,1,k) + (Ey(i,1,k+1) - Ey(i,1,k)) * dt_invdz
		end do
		do i = 1, nx - 1
			By(i,1,k) = By(i,1,k) - (Ex(i,1,k+1) - Ex(i,1,k)) * dt_invdz &
									+ (Ez(i+1,1,k) - Ez(i,1,k)) * dt_invdx
			Bz(i,1,k) = Bz(i,1,k) - (Ey(i+1,1,k) - Ey(i,1,k)) * dt_invdx
		end do
	end do
	call Ghost_Cells(Bx, 'Bx', domain); call Ghost_Cells(By, 'By', domain); call Ghost_Cells(Bz, 'Bz', domain)
end subroutine Update_B_FDTD
!*********************************************************************
subroutine Update_E_FDTD(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, domain)	! 电场更新器
    use physicals_parameter; use Ghost_Cells_Module
	use mesh_parameter, only: nx, nz, inv_dx, inv_dz
	use Domain_Decomposition_Module, only: domain_decomposition
    implicit none
	real(8), intent(inout) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
    real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:)
    real(8), intent(in) :: Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    integer :: i, k, k_backward_start, k_interior_end, k_forward_end
    real(8) :: curl_invdx, curl_invdz

	curl_invdx = coef_curl * inv_dx; curl_invdz = coef_curl * inv_dz
    k_backward_start = 1
    if (domain%k_start == 1) k_backward_start = 2
    k_interior_end = domain%nz_local
    if (domain%k_end == nz) k_interior_end = domain%nz_local - 1
    k_forward_end = domain%nz_local
    if (domain%k_end == nz) k_forward_end = domain%nz_local - 1

    ! Ampere 方程在 XZ 平面上的形式：
    ! dEx/dt=-c2*dBy/dz-Jx/eps0,
    ! dEy/dt=c2*(dBx/dz-dBz/dx)-Jy/eps0,
    ! dEz/dt=-c2*dBy/dx-Jz/eps0。
    do k = k_backward_start, k_interior_end; do i = 1, nx - 1
        Ex(i,1,k) = Ex(i,1,k) - (By(i,1,k) - By(i,1,k-1))*curl_invdz - Jx(i,1,k)*coef_current
    end do; end do
    do k = k_backward_start, k_interior_end; do i = 2, nx - 1
        Ey(i,1,k) = Ey(i,1,k) + (Bx(i,1,k) - Bx(i,1,k-1))*curl_invdz &
                    - (Bz(i,1,k) - Bz(i-1,1,k))*curl_invdx - Jy(i,1,k)*coef_current
    end do; end do
    do k = 1, k_forward_end; do i = 2, nx - 1
        Ez(i,1,k) = Ez(i,1,k) - (By(i,1,k) - By(i-1,1,k))*curl_invdx - Jz(i,1,k)*coef_current
    end do; end do

    Ey(1, :, :) = 0.0d0; Ey(nx, :, :) = 0.0d0
    Ez(1, :, :) = 0.0d0; Ez(nx, :, :) = 0.0d0
    call Ghost_Cells(Ex, 'Ex', domain); call Ghost_Cells(Ey, 'Ey', domain); call Ghost_Cells(Ez, 'Ez', domain)
end subroutine Update_E_FDTD
!*********************************************************************
end program pic_equilibrium
