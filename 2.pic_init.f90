!***************************************************************
!	Coded By Siyu Zhang     May. 19, 2026
!
!   Full-PIC 平衡采样 高温等离子体物理学课题用
!
!	基于 Full-PIC-init_v3.0 版本修改 (Siyu Zhang, November. 16, 2025)
!*********************************************************************
module mesh_parameter
implicit none
	integer, parameter :: nx = 257, ny = 257, nz = 3		! 网格数
    real(8), parameter :: x_min = -6.4d-2		! 网格的物理边界
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
end module physicals_parameter
!*********************************************************************
module pic_types
use physicals_parameter, only: q, pi
implicit none
	type :: particle
		real(8) :: x, y, z		! 粒子位置 (m)
		real(8) :: vx, vy, vz	! 粒子速度 (m/s)
		real(8) :: charge, mass	! 粒子电荷 (C) 和质量 (kg)
		real(8) :: weight		! 超粒子权重
	end type particle
contains
!----------------------------------
function sample_maxwell(temp_ev, mass) result(v)
implicit none
	real(8), intent(in) :: temp_ev, mass
	real(8) :: sigma, u1, u2, u3, u4
	real(8) :: v(3)
	real(8) :: term1, term2
	real(8) :: temp_joule
	
	temp_joule = temp_ev * q
    sigma = sqrt(temp_joule / mass)
	
	! 生成第一对正态分布随机数
	call random_number(u1)
	call random_number(u2)
	term1 = sigma * sqrt(-2d0 * log(u1))
	v(1) = term1 * cos(2d0 * pi * u2)
	v(2) = term1 * sin(2d0 * pi * u2)
	
	! 生成第二对
	call random_number(u3)
	call random_number(u4)
	term2 = sigma * sqrt(-2d0 * log(u3))
	v(3) = term2 * cos(2d0 * pi * u4)
end function sample_maxwell
end module pic_types
!*********************************************************************
program pic_equilibrium
use mpi; use mesh_parameter; use physicals_parameter; use pic_types
implicit none
	integer :: rank, nprocs, ierr
    real(8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:), &
							Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), &
							dens(:,:,:), prs(:,:,:), psi(:,:,:)
	integer :: j, p, idyn
	real(8) :: w_sum_ion = 0d0, w_sum_elec = 0d0	! 权重和，电中性检验用
	integer :: nIon  = 50000000		! 超离子数
	integer :: nElec = 50000000		! 超电子数
	integer :: nIon_local, nElec_local		! 各进程分配粒子用
	integer :: ion_offset, elec_offset
	integer :: base_Ion, extra_Ion, base_Elec, extra_Elec
	type(particle), allocatable :: ion_local(:), elec_local(:)
	integer :: seed_size, i			! 随机数用
	integer, allocatable :: seed(:)
	integer :: base_seed = 13579    ! 种子

	! 初始化 MPI
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
	
	! 确定性随机数种子数组生成
	call random_seed(size = seed_size)
	allocate(seed(seed_size))
	do i = 1, seed_size
		seed(i) = int(mod(dble(base_seed) + 1299721.0d0 * rank + 97.0d0 * i, 2147483647.0d0))
	end do
	call random_seed(put = seed)
	deallocate(seed)

	if (rank == 0) print *, '----------------- Start initialization -----------------'
	
	! 计算本地粒子数和偏移量
	base_Ion = nIon / nprocs
	extra_Ion = mod(nIon, nprocs)
	if (rank < extra_Ion) then
		nIon_local = base_Ion + 1
		ion_offset = rank * (base_Ion + 1)
	else
		nIon_local = base_Ion
		ion_offset = extra_Ion * (base_Ion + 1) + (rank - extra_Ion) * base_Ion
	end if

	base_Elec = nElec / nprocs
	extra_Elec = mod(nElec, nprocs)
	if (rank < extra_Elec) then
		nElec_local = base_Elec + 1
		elec_offset = rank * (base_Elec + 1)
	else
		nElec_local = base_Elec
		elec_offset = extra_Elec * (base_Elec + 1) + (rank - extra_Elec) * base_Elec
	end if
	
	! 进程数检查
	if (rank==0) print *, '[Process count] NCPU =', nprocs
	
	! 分配场数组,粒子数组
	allocate(Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz), &
			Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz), &
			dens(nx,ny,nz), prs(nx,ny,nz), psi(nx,ny,nz))
	allocate(ion_local(nIon_local), elec_local(nElec_local))

	! 加载初始MHD平衡场(t=0时刻的场),并广播
	if (rank == 0) call load_equilibrium(Bx, By, Bz, dens, prs, Ex, Ey, Ez, psi)
	call MPI_Bcast(Bx, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(By, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Bz, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ex, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ey, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ez, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dens, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(prs, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(psi, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)	

	! 初始化离子
	call init_ions(ion_local, nIon_local, nIon, Bx, By, Bz, dens, prs, w_sum_ion, rank, ierr)
	
	! 随机数种子加偏移
	call random_seed(size = seed_size)
	allocate(seed(seed_size))
	do i = 1, seed_size
		seed(i) = int(mod(dble(base_seed + 99999) + 1299721.0d0 * rank + 97.0d0 * i, 2147483647.0d0))
	end do
	call random_seed(put = seed)
	deallocate(seed)
	
	! 初始化电子
	call init_electrons(elec_local, nElec_local, nElec, Bx, By, Bz, Ex, Ey, Ez, dens, prs, w_sum_elec, rank, ierr)
	
	! 电中性检查
	if (rank==0) print *, '[charge neutrality] Q_total =', q*(w_sum_ion - w_sum_elec)

	! 粒子数据写入文件
	if (rank == 0) print *, '----------------- Start writing particle data -----------------'
	call write_particle_data(ion_local, elec_local, nIon_local, nElec_local, rank, nprocs, ierr)
	call MPI_Barrier(MPI_COMM_WORLD, ierr)

	! 粒子位置写入bin文件
	if (rank == 0) print *, '----------------- Start writing position data -----------------'
	call dump_dynamic_binary(ion_local, nIon_local, "ion", 0)
	call dump_dynamic_binary(elec_local, nElec_local, "elec", 0)
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (rank == 0) print *, '----------------- End -------------------'
	call MPI_Finalize(ierr)
contains
!*********************************************************************
subroutine load_equilibrium(Bx, By, Bz, dens, prs, Ex, Ey, Ez, psi)		! 读背景场文件
use mesh_parameter, only: nx, ny, nz
implicit none
    real(8), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz), &
							Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz), &
							dens(nx,ny,nz), prs(nx,ny,nz), psi(nx,ny,nz)
    integer :: io_stat
    open(unit=10,file='out_3d/p_3d.bin',form='unformatted',access='stream',status='old',iostat=io_stat)
    read(10) prs; close(10)   
    open(unit=10, file='out_3d/psi_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) psi; close(10)
    open(unit=10, file='out_3d/den_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) dens; close(10)
    open(unit=10, file='out_3d/Bx_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) Bx; close(10)
    open(unit=10, file='out_3d/By_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) By; close(10)
    open(unit=10, file='out_3d/Bz_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) Bz; close(10)
    open(unit=10, file='out_3d/Ex_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) Ex; close(10)
    open(unit=10, file='out_3d/Ey_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) Ey; close(10)
    open(unit=10, file='out_3d/Ez_3d.bin', form='unformatted', access='stream', status='old', iostat=io_stat)
    read(10) Ez; close(10)
end subroutine load_equilibrium
!*********************************************************************
subroutine init_ions(ion_local, nI_local, nI_global, Bx, By, Bz, dens, prs, w_sum_ion_global, rank, ierr)	! 离子初始化
use mesh_parameter; use physicals_parameter; use pic_types
implicit none
	integer, intent(in) :: nI_local, nI_global, rank
	type(particle), intent(out) :: ion_local(nI_local)
	real(8), intent(in) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz), dens(nx,ny,nz), prs(nx,ny,nz)
	real(8), intent(out) :: w_sum_ion_global	! 总权重
	integer :: i, ix, iy, iz
	real(8) :: xx, yy, zz, u	! 随机数 (0,1)
	real(8) :: x_phys, y_phys, z_phys
	real(8) :: v3(3), dens_max
	real(8) :: vol_cell, n_real_i, weight_const, weighting_ion
	real(8), allocatable :: n_chk(:,:,:), n_chk_cell(:,:,:) ! 3D 沉积检查数组
	real(8) :: err_threshold, dens_floor, err_max = 0		! 沉积检查用
	real(8) :: err_sum_sq = 0.0d0   	! 误差平方和
	real(8) :: rmse_val = 0.0d0     	! 最终的RMSE值
	real(8) :: w_sum_ion_local = 0.0d0	! 本地权重和
	real(8) :: mass = 25.0d0 * 9.1093837139d-31, temp = 1000.0d0	! 缩放物理参数：质量比25，温度1000 eV
	real(8) :: alpha, beta, gamma, dens_interpolated	! 三线性插值用
	real(8) :: x_grid, y_grid, z_grid	! CIC用
	real(8) :: Bx_interp, By_interp, Bz_interp, Ex_interp, Ey_interp, Ez_interp
    real(8) :: B_sq, B_mag, v_th, drift_mag, limit_factor
    real(8) :: grad_Px, grad_Py, grad_Pz, Pi_frac
    real(8) :: vE_x, vE_y, vE_z, vDia_x, vDia_y, vDia_z
    integer :: ix_m, ix_p, iy_m, iy_p, iz_m, iz_p
	integer :: n_cells_checked = 0		! 被计入统计的网格单元数
	integer :: ix_idx, iy_idx, iz_idx	! 粒子所在网格索引
	integer, intent(inout) :: ierr

	! 计算真实粒子数
	n_real_i = 0d0
	if (rank == 0) then
        do ix = 1, nx-1
            do iy = 1, ny-1
                do iz = 1, nz-1 
                    vol_cell = dx * dy * dz		! 笛卡尔单元体积 
                    n_real_i = n_real_i + dens(ix,iy,iz) * vol_cell		! 使用单元中心或左下角密度
                end do
            end do
        end do
        print *, '[Background field based] ion =', n_real_i
    end if
    call MPI_Bcast(n_real_i, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

	weight_const = n_real_i / dble(nI_global)	! 常数权重
	dens_max = maxval(dens)		! 预备拒绝采样参数和沉积检查用

	do i = 1, nI_local
		do	! 拒绝采样+CIC
			call random_number(xx)
			call random_number(yy)
			call random_number(zz)
			call random_number(u)
			x_phys = x_min + xx * (dble(nx-1) * dx)
			y_phys = y_min + yy * (dble(ny-1) * dy)
			z_phys = z_min + zz * (dble(nz-1) * dz)			
			
			! 计算粒子所在单元的左下角索引
			ix_idx = int((x_phys - x_min) / dx) + 1
			iy_idx = int((y_phys - y_min) / dy) + 1
			iz_idx = int((z_phys - z_min) / dz) + 1
			
			! 避免越界
			if (ix_idx<1 .or. ix_idx>=nx .or. iy_idx<1 .or. iy_idx>=ny .or. iz_idx<1 .or. iz_idx>=nz) cycle

			! 计算在单元内的归一化相对位置
			x_grid = dble(ix_idx-1) * dx + x_min
			y_grid = dble(iy_idx-1) * dy + y_min
			z_grid = dble(iz_idx-1) * dz + z_min
            
			alpha = (x_phys - x_grid) / dx
			beta  = (y_phys - y_grid) / dy
            gamma = (z_phys - z_grid) / dz

			alpha = max(0.0d0, min(1.0d0, alpha))
			beta  = max(0.0d0, min(1.0d0, beta))
            gamma = max(0.0d0, min(1.0d0, gamma))

			! 对背景密度场三线性插值
            dens_interpolated = dens(ix_idx, iy_idx, iz_idx) * (1.d0-alpha) * (1.d0-beta) * (1.d0-gamma) + &
								dens(ix_idx+1, iy_idx, iz_idx) * alpha * (1.d0-beta) * (1.d0-gamma) + &
								dens(ix_idx, iy_idx+1, iz_idx) * (1.d0-alpha) * beta * (1.d0-gamma) + &
								dens(ix_idx, iy_idx, iz_idx+1) * (1.d0-alpha) * (1.d0-beta) * gamma + &
								dens(ix_idx+1, iy_idx+1, iz_idx) * alpha * beta * (1.d0-gamma) + &
								dens(ix_idx+1, iy_idx, iz_idx+1) * alpha * (1.d0-beta) * gamma + &
								dens(ix_idx, iy_idx+1, iz_idx+1) * (1.d0-alpha) * beta * gamma + &
								dens(ix_idx+1, iy_idx+1, iz_idx+1) * alpha * beta * gamma
			if (u < dens_interpolated / dens_max) exit
		end do
		
		! 存储坐标
		ion_local(i)%x = x_phys
		ion_local(i)%y = y_phys
		ion_local(i)%z = z_phys
		v3 = sample_maxwell(temp, mass)

		! 场插值
		Bx_interp = trilinear_interp(Bx, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		By_interp = trilinear_interp(By, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Bz_interp = trilinear_interp(Bz, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Ex_interp = trilinear_interp(Ex, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Ey_interp = trilinear_interp(Ey, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Ez_interp = trilinear_interp(Ez, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		
		B_sq = Bx_interp**2 + By_interp**2 + Bz_interp**2 + 1.0d-20
		B_mag = sqrt(B_sq)

		! 计算离子压力梯度 (假设 Ti = Te，离子压力占总压力的 0.5)
		Pi_frac = 0.5d0
		ix_m = max(1, ix_idx-1); ix_p = min(nx, ix_idx+1)
		iy_m = max(1, iy_idx-1); iy_p = min(ny, iy_idx+1)
		iz_m = max(1, iz_idx-1); iz_p = min(nz, iz_idx+1)
		
		grad_Px = Pi_frac * (prs(ix_p, iy_idx, iz_idx) - prs(ix_m, iy_idx, iz_idx)) / (dble(ix_p - ix_m) * dx)
		grad_Py = Pi_frac * (prs(ix_idx, iy_p, iz_idx) - prs(ix_idx, iy_m, iz_idx)) / (dble(iy_p - iy_m) * dy)
		grad_Pz = Pi_frac * (prs(ix_idx, iy_idx, iz_p) - prs(ix_idx, iy_idx, iz_m)) / (dble(iz_p - iz_m) * dz)

		! 计算漂移速度 E x B 漂移
		vE_x = (Ey_interp * Bz_interp - Ez_interp * By_interp) / B_sq
		vE_y = (Ez_interp * Bx_interp - Ex_interp * Bz_interp) / B_sq
		vE_z = (Ex_interp * By_interp - Ey_interp * Bx_interp) / B_sq
		
		! 离子抗磁漂移 (B x grad(P) / (q n B^2))
		vDia_x = (By_interp * grad_Pz - Bz_interp * grad_Py) / (q * dens_interpolated * B_sq)
		vDia_y = (Bz_interp * grad_Px - Bx_interp * grad_Pz) / (q * dens_interpolated * B_sq)
		vDia_z = (Bx_interp * grad_Py - By_interp * grad_Px) / (q * dens_interpolated * B_sq)

		! 漂移速度限幅 (防止磁零点奇异性导致程序崩溃)
		drift_mag = sqrt((vE_x + vDia_x)**2 + (vE_y + vDia_y)**2 + (vE_z + vDia_z)**2)
		v_th = sqrt(2.0d0 * temp * q / mass)
		limit_factor = 1.0d0
		if (drift_mag > 2.0d0 * v_th) then
			limit_factor = (2.0d0 * v_th) / drift_mag
		end if

		! 组合最终速度：热速度 + 限制后的漂移速度
		ion_local(i)%vx = v3(1) + (vE_x + vDia_x) * limit_factor
		ion_local(i)%vy = v3(2) + (vE_y + vDia_y) * limit_factor
		ion_local(i)%vz = v3(3) + (vE_z + vDia_z) * limit_factor

		! 物性与权重
		ion_local(i)%mass = mass
		ion_local(i)%charge = q
		ion_local(i)%weight = weight_const
	end do
	
	! 离子权重
	w_sum_ion_global = 0.0d0
	if (nI_local > 0) then
		weighting_ion = ion_local(1)%weight
	else
		weighting_ion = 0
	end if
	if (rank == 0) print *, '[Hyperparticle Weighting] ion =', weighting_ion
	
	do i = 1, nI_local
		w_sum_ion_local = w_sum_ion_local + ion_local(i)%weight
	end do
    call MPI_Reduce(w_sum_ion_local, w_sum_ion_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

	! 密度回算,粒子沉积检查
	allocate(n_chk(nx,ny,nz)); allocate(n_chk_cell(nx,ny,nz))
	n_chk = 0d0; n_chk_cell = 0d0
	
	call deposit_weights(ion_local, nI_local, n_chk)	! 各进程粒子权重沉积到本地数组
	call MPI_Allreduce(MPI_IN_PLACE, n_chk, nx*ny*nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

	! 节点到单元中心插值
	if (rank == 0) then
		vol_cell = dx * dy * dz
		if (vol_cell <= 1.0d-30) then
            print *, "FATAL ERROR in init_ions: Cell volume is zero or negative."
            call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        end if
		
		do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    ! 单元中心 (8个节点平均)
                    n_chk_cell(ix, iy, iz) = 0.125d0 * (n_chk(ix, iy, iz) + n_chk(ix+1, iy, iz) + &
														n_chk(ix, iy+1, iz) + n_chk(ix, iy, iz+1) + &
														n_chk(ix+1, iy+1, iz) + n_chk(ix+1, iy, iz+1) + &
														n_chk(ix, iy+1, iz+1) + n_chk(ix+1, iy+1, iz+1))
					! 单元中心权重除以单元体积，得到单元中心密度
                    n_chk_cell(ix, iy, iz) = n_chk_cell(ix, iy, iz) / vol_cell
				end do
            end do
		end do
		dens_floor = 1.0d-1 * dens_max	! 密度差阈值	
		err_max = 0.0d0
		err_sum_sq = 0.0d0
		n_cells_checked = 0
		do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    if (dens(ix,iy,iz) > dens_floor) then
                        err_threshold = abs(n_chk_cell(ix,iy,iz) - dens(ix,iy,iz)) / dens(ix,iy,iz)
                        if (err_threshold > err_max) err_max = err_threshold
                        
                        err_sum_sq = err_sum_sq + (n_chk_cell(ix,iy,iz) - dens(ix,iy,iz))**2
                        n_cells_checked = n_cells_checked + 1
                    end if
                end do
            end do
		end do	
		rmse_val = 0.0d0
		if (n_cells_checked > 0) rmse_val = sqrt(err_sum_sq / dble(n_cells_checked))
		print *, '[ION] init inaccuracies =', err_max
		print *, '[ION] RMSE of density =', rmse_val
	end if
	deallocate(n_chk); deallocate(n_chk_cell)
end subroutine init_ions
!*********************************************************************
subroutine init_electrons(elec_local, nE_local, nE_global, Bx, By, Bz, Ex, Ey, Ez, dens, prs, w_sum_elec_global, rank, ierr)	! 电子初始化
use mesh_parameter; use physicals_parameter; use pic_types
implicit none
	integer, intent(in) :: nE_local, nE_global, rank
	type(particle), intent(out) :: elec_local(nE_local)
	real(8), intent(in) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz), &
						   Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz), &
						   dens(nx,ny,nz), prs(nx,ny,nz)
	real(8), intent(out) :: w_sum_elec_global	
	integer, intent(inout) :: ierr
	integer :: i, ix, iy, iz
	real(8) :: xx, yy, zz, u
	real(8) :: x_phys, y_phys, z_phys
	real(8) :: v3(3), dens_max
	real(8) :: vol_cell, n_real_e, weight_const, weighting_elec
	real(8), allocatable :: n_chk(:,:,:), n_chk_cell(:,:,:)
	real(8) :: err_threshold, dens_floor, err_max = 0
	real(8) :: err_sum_sq = 0.0d0
	real(8) :: rmse_val = 0.0d0
	real(8) :: w_sum_elec_local = 0.0d0
	real(8) :: alpha, beta, gamma, dens_interpolated
    real(8) :: x_grid, y_grid, z_grid
    real(8) :: Bx_interp, By_interp, Bz_interp, Ex_interp, Ey_interp, Ez_interp	! 存储插值后的场分量
	real(8) :: B_sq, drift_x, drift_y, drift_z
	real(8) :: B_mag, v_th, drift_mag, limit_factor
    real(8) :: grad_Px, grad_Py, grad_Pz, Pe_frac
    real(8) :: vE_x, vE_y, vE_z, vDia_x, vDia_y, vDia_z
    integer :: ix_m, ix_p, iy_m, iy_p, iz_m, iz_p
    real(8) :: mass = 9.1093837139d-31, temp = 1000.0d0
    integer :: ix_idx, iy_idx, iz_idx
	integer :: n_cells_checked = 0

	n_real_e = 0d0
	if (rank == 0) then
		do ix = 1, nx-1
            do iy = 1, ny-1
                do iz = 1, nz-1
                    vol_cell = dx * dy * dz
                    n_real_e = n_real_e + dens(ix,iy,iz) * vol_cell
                end do
            end do
        end do
		print *, '[Background field based] elec =', n_real_e
	end if
	call MPI_Bcast(n_real_e, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

	weight_const = n_real_e / dble(nE_global)
	dens_max = maxval(dens)

	do i = 1, nE_local
		do
			call random_number(xx)
			call random_number(yy)
			call random_number(zz)
			call random_number(u)
			
			x_phys = x_min + xx * (dble(nx-1) * dx)
			y_phys = y_min + yy * (dble(ny-1) * dy)
			z_phys = z_min + zz * (dble(nz-1) * dz)
			
			ix_idx = int((x_phys - x_min) / dx) + 1
			iy_idx = int((y_phys - y_min) / dy) + 1
			iz_idx = int((z_phys - z_min) / dz) + 1
			
			if (ix_idx<1 .or. ix_idx>=nx .or. iy_idx<1 .or. iy_idx>=ny .or. iz_idx<1 .or. iz_idx>=nz) cycle
			
			x_grid = dble(ix_idx-1) * dx + x_min
			y_grid = dble(iy_idx-1) * dy + y_min
			z_grid = dble(iz_idx-1) * dz + z_min
			
			alpha = (x_phys - x_grid) / dx
			beta = (y_phys - y_grid) / dy
            gamma = (z_phys - z_grid) / dz
			
			alpha = max(0.0d0, min(1.0d0, alpha))
			beta = max(0.0d0, min(1.0d0, beta))
            gamma = max(0.0d0, min(1.0d0, gamma))
			
			dens_interpolated = dens(ix_idx, iy_idx, iz_idx) * (1.d0-alpha) * (1.d0-beta) * (1.d0-gamma) + &
								dens(ix_idx+1, iy_idx, iz_idx) * alpha * (1.d0-beta) * (1.d0-gamma) + &
								dens(ix_idx, iy_idx+1, iz_idx) * (1.d0-alpha) * beta * (1.d0-gamma) + &
								dens(ix_idx, iy_idx, iz_idx+1) * (1.d0-alpha) * (1.d0-beta) * gamma + &
								dens(ix_idx+1, iy_idx+1, iz_idx) * alpha * beta * (1.d0-gamma) + &
								dens(ix_idx+1, iy_idx, iz_idx+1) * alpha * (1.d0-beta) * gamma + &
								dens(ix_idx, iy_idx+1, iz_idx+1) * (1.d0-alpha) * beta * gamma + &
								dens(ix_idx+1, iy_idx+1, iz_idx+1) * alpha * beta * gamma
				
			if (u < dens_interpolated / dens_max) exit
		end do
		
		elec_local(i)%x = x_phys
        elec_local(i)%y = y_phys
        elec_local(i)%z = z_phys	
		v3 = sample_maxwell(temp, mass)

		! 场插值
		Bx_interp = trilinear_interp(Bx, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		By_interp = trilinear_interp(By, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Bz_interp = trilinear_interp(Bz, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Ex_interp = trilinear_interp(Ex, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Ey_interp = trilinear_interp(Ey, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		Ez_interp = trilinear_interp(Ez, ix_idx, iy_idx, iz_idx, alpha, beta, gamma)
		
		B_sq = Bx_interp**2 + By_interp**2 + Bz_interp**2 + 1.0d-20
		
		! 电子压力分量
		Pe_frac = 0.5d0
		ix_m = max(1, ix_idx-1); ix_p = min(nx, ix_idx+1)
		iy_m = max(1, iy_idx-1); iy_p = min(ny, iy_idx+1)
		iz_m = max(1, iz_idx-1); iz_p = min(nz, iz_idx+1)
		
		grad_Px = Pe_frac * (prs(ix_p, iy_idx, iz_idx) - prs(ix_m, iy_idx, iz_idx)) / (dble(ix_p - ix_m) * dx)
		grad_Py = Pe_frac * (prs(ix_idx, iy_p, iz_idx) - prs(ix_idx, iy_m, iz_idx)) / (dble(iy_p - iy_m) * dy)
		grad_Pz = Pe_frac * (prs(ix_idx, iy_idx, iz_p) - prs(ix_idx, iy_idx, iz_m)) / (dble(iz_p - iz_m) * dz)

		! E x B 漂移
		vE_x = (Ey_interp * Bz_interp - Ez_interp * By_interp) / B_sq
		vE_y = (Ez_interp * Bx_interp - Ex_interp * Bz_interp) / B_sq
		vE_z = (Ex_interp * By_interp - Ey_interp * Bx_interp) / B_sq
		
		! 电子抗磁漂移
		vDia_x = -(By_interp * grad_Pz - Bz_interp * grad_Py) / (q * dens_interpolated * B_sq)
		vDia_y = -(Bz_interp * grad_Px - Bx_interp * grad_Pz) / (q * dens_interpolated * B_sq)
		vDia_z = -(Bx_interp * grad_Py - By_interp * grad_Px) / (q * dens_interpolated * B_sq)

		! 限幅
		drift_mag = sqrt((vE_x + vDia_x)**2 + (vE_y + vDia_y)**2 + (vE_z + vDia_z)**2)
		v_th = sqrt(2.0d0 * temp * q / mass)
		limit_factor = 1.0d0
		if (drift_mag > 2.0d0 * v_th) then
			limit_factor = (2.0d0 * v_th) / drift_mag
		end if

		elec_local(i)%vx = v3(1) + (vE_x + vDia_x) * limit_factor
		elec_local(i)%vy = v3(2) + (vE_y + vDia_y) * limit_factor
		elec_local(i)%vz = v3(3) + (vE_z + vDia_z) * limit_factor
		elec_local(i)%mass = mass
		elec_local(i)%charge = -q
		elec_local(i)%weight = weight_const
	end do
	
	if (nE_local > 0) then
		weighting_elec = elec_local(1)%weight
	else
		weighting_elec = 0
	end if
    if (rank == 0) print *, '[Hyperparticle Weighting] ele =', weighting_elec

	w_sum_elec_global = 0.0d0
    do i = 1, nE_local
        w_sum_elec_local = w_sum_elec_local + elec_local(i)%weight
    end do
    call MPI_Reduce(w_sum_elec_local, w_sum_elec_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

	allocate(n_chk(nx,ny,nz)); allocate(n_chk_cell(nx,ny,nz))
	n_chk = 0d0; n_chk_cell = 0d0
	
    call deposit_weights(elec_local, nE_local, n_chk)
	call MPI_Allreduce(MPI_IN_PLACE, n_chk, nx*ny*nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

	if (rank == 0) then
		vol_cell = dx * dy * dz
        if (vol_cell <= 1.0d-30) then
            print *, "FATAL ERROR in init_electrons: Cell volume is zero or negative."
            call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        end if
		
		do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    n_chk_cell(ix, iy, iz) = 0.125d0 * (n_chk(ix, iy, iz) + n_chk(ix+1, iy, iz) + &
														n_chk(ix, iy+1, iz) + n_chk(ix, iy, iz+1) + &
														n_chk(ix+1, iy+1, iz) + n_chk(ix+1, iy, iz+1) + &
														n_chk(ix, iy+1, iz+1) + n_chk(ix+1, iy+1, iz+1))
					n_chk_cell(ix, iy, iz) = n_chk_cell(ix, iy, iz) / vol_cell
				end do
            end do
		end do
		dens_floor = 1.0d-1 * dens_max
		err_max = 0.0d0
		err_sum_sq = 0.0d0
		n_cells_checked = 0
		do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    if (dens(ix,iy,iz) > dens_floor) then
                        err_threshold = abs(n_chk_cell(ix,iy,iz) - dens(ix,iy,iz)) / dens(ix,iy,iz)
                        if (err_threshold > err_max) err_max = err_threshold
                        
                        err_sum_sq = err_sum_sq + (n_chk_cell(ix,iy,iz) - dens(ix,iy,iz))**2
                        n_cells_checked = n_cells_checked + 1
                    end if
                end do
            end do
		end do

        rmse_val = 0.0d0
        if (n_cells_checked > 0) rmse_val = sqrt(err_sum_sq / dble(n_cells_checked))
        print *, '[elec] init inaccuracies =', err_max
        print *, '[elec] RMSE of density =', rmse_val
    end if
	deallocate(n_chk); deallocate(n_chk_cell)
end subroutine init_electrons
!*********************************************************************
function trilinear_interp(field, ix, iy, iz, a, b, g) result(field_interp)	! BxE用场插值
	implicit none
	real(8), intent(in) :: field(nx,ny,nz)
	integer, intent(in) :: ix, iy, iz
	real(8), intent(in) :: a, b, g
	real(8) :: field_interp
	field_interp = field(ix, iy, iz) * (1.d0-a) * (1.d0-b) * (1.d0-g) + &
				   field(ix+1, iy, iz) * a * (1.d0-b) * (1.d0-g) + &
				   field(ix, iy+1, iz) * (1.d0-a) * b * (1.d0-g) + &
				   field(ix, iy, iz+1) * (1.d0-a) * (1.d0-b) * g + &
				   field(ix+1, iy+1, iz) * a * b * (1.d0-g) + &
				   field(ix+1, iy, iz+1) * a * (1.d0-b) * g + &
				   field(ix, iy+1, iz+1) * (1.d0-a) * b * g + &
				   field(ix+1, iy+1, iz+1) * a * b * g
end function trilinear_interp
!*********************************************************************
subroutine deposit_weights(p, npart, w_grid)	! 粒子权重沉积到网格节点上
use mesh_parameter; use pic_types
implicit none
    type(particle), intent(in) :: p(:)
    integer, intent(in) :: npart
    real(8), intent(inout) :: w_grid(nx, ny, nz)
    integer :: k, ix, iy, iz
    real(8) :: x_phys, y_phys, z_phys
	real(8) :: alpha, beta, gamma		! 3D 相对位置
    real(8) :: w000, w100, w010, w001, w110, w101, w011, w111	! 8 个角的权重
    real(8) :: x_i0, y_j0, z_k0			! 单元左下角坐标

    ! 沉积加权粒子数到网格节点
    do k = 1, npart
		x_phys = p(k)%x
        y_phys = p(k)%y
        z_phys = p(k)%z

		! 计算粒子所在单元的左下角索引
        ix = int((x_phys - x_min) / dx) + 1
        iy = int((y_phys - y_min) / dy) + 1
        iz = int((z_phys - z_min) / dz) + 1

		if (ix<1 .or. ix>nx-1 .or. iy<1 .or. iy>ny-1 .or. iz<1 .or. iz>nz-1) cycle

		! 获取单元左下角的坐标
        x_i0 = dble(ix-1) * dx + x_min
        y_j0 = dble(iy-1) * dy + y_min
        z_k0 = dble(iz-1) * dz + z_min
		
		! 计算归一化相对位置
        alpha = (x_phys - x_i0) / dx
        beta = (y_phys - y_j0) / dy
        gamma = (z_phys - z_k0) / dz
        alpha = max(0.0d0, min(1.0d0, alpha))
        beta = max(0.0d0, min(1.0d0, beta))
        gamma = max(0.0d0, min(1.0d0, gamma))	
		
		! 计算三线性插值的权重
        w000 = (1.0d0 - alpha) * (1.0d0 - beta) * (1.0d0 - gamma)
        w100 = alpha * (1.0d0 - beta) * (1.0d0 - gamma)
        w010 = (1.0d0 - alpha) * beta * (1.0d0 - gamma)
        w001 = (1.0d0 - alpha) * (1.0d0 - beta) * gamma
        w110 = alpha * beta * (1.0d0 - gamma)
        w101 = alpha * (1.0d0 - beta) * gamma
        w011 = (1.0d0 - alpha) * beta * gamma
        w111 = alpha * beta * gamma

        ! 沉积到8个节点
        w_grid(ix, iy, iz) = w_grid(ix, iy, iz) + p(k)%weight * w000
        w_grid(ix+1, iy, iz) = w_grid(ix+1, iy, iz) + p(k)%weight * w100
        w_grid(ix, iy+1, iz) = w_grid(ix, iy+1, iz) + p(k)%weight * w010
        w_grid(ix, iy, iz+1) = w_grid(ix, iy, iz+1) + p(k)%weight * w001
        w_grid(ix+1, iy+1, iz) = w_grid(ix+1, iy+1, iz) + p(k)%weight * w110
        w_grid(ix+1, iy, iz+1) = w_grid(ix+1, iy, iz+1) + p(k)%weight * w101
        w_grid(ix, iy+1, iz+1) = w_grid(ix, iy+1, iz+1) + p(k)%weight * w011
        w_grid(ix+1, iy+1, iz+1) = w_grid(ix+1, iy+1, iz+1) + p(k)%weight * w111	
    end do
end subroutine deposit_weights
!*********************************************************************
subroutine write_particle_data(ion_local, elec_local, nIon_local, nElec_local, rank, nprocs, ierr)	! 粒子初始化信息写入文件
use pic_types, only: particle; use mpi
implicit none
    type(particle), intent(in) :: ion_local(:), elec_local(:)
    integer, intent(in) :: nIon_local, nElec_local, rank, nprocs
	integer, intent(inout) :: ierr
    integer :: p
    character(len=64) :: filename
	character(len=10) :: open_status
	
    ! 各进程轮流写入
    do p = 0, nprocs - 1
        if (rank == p) then
            if (p == 0) then		! rank 0 替换，其他追加
                open_status = 'replace'
            else
                open_status = 'append'
            end if

			! 因为用于启动PIC，所以不做权重检查，写入全部粒子(SoA 布局)
			write(filename, '(A,A)') 'init_particles/ion_x', '.bin'		! 离子 x 位置
            call write_component_binary(filename, ion_local, nIon_local, 'x', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/ion_y', '.bin'		! 离子 y 位置
            call write_component_binary(filename, ion_local, nIon_local, 'y', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/ion_z', '.bin'  	! 离子 z 位置
            call write_component_binary(filename, ion_local, nIon_local, 'z', open_status, rank, MPI_COMM_WORLD, ierr)    
            write(filename, '(A,A)') 'init_particles/ion_vx', '.bin' 	! 离子 x 速度
            call write_component_binary(filename, ion_local, nIon_local, 'vx', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/ion_vy', '.bin' 	! 离子 y 速度
            call write_component_binary(filename, ion_local, nIon_local, 'vy', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/ion_vz', '.bin' 	! 离子 z 速度
            call write_component_binary(filename, ion_local, nIon_local, 'vz', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/ion_weight', '.bin'	! 离子权重
            call write_component_binary(filename, ion_local, nIon_local, 'weight', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/elec_x', '.bin' 	! 电子 x 位置
            call write_component_binary(filename, elec_local, nElec_local, 'x', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/elec_y', '.bin' 	! 电子 y 位置
            call write_component_binary(filename, elec_local, nElec_local, 'y', open_status, rank, MPI_COMM_WORLD, ierr)
            write(filename, '(A,A)') 'init_particles/elec_z', '.bin' 	! 电子 z 位置
            call write_component_binary(filename, elec_local, nElec_local, 'z', open_status, rank, MPI_COMM_WORLD, ierr)           
            write(filename, '(A,A)') 'init_particles/elec_vx', '.bin' 	! 电子 x 速度
            call write_component_binary(filename, elec_local, nElec_local, 'vx', open_status, rank, MPI_COMM_WORLD, ierr)            
            write(filename, '(A,A)') 'init_particles/elec_vy', '.bin' 	! 电子 y 速度
            call write_component_binary(filename, elec_local, nElec_local, 'vy', open_status, rank, MPI_COMM_WORLD, ierr)            
            write(filename, '(A,A)') 'init_particles/elec_vz', '.bin' 	! 电子 z 速度
            call write_component_binary(filename, elec_local, nElec_local, 'vz', open_status, rank, MPI_COMM_WORLD, ierr)            
            write(filename, '(A,A)') 'init_particles/elec_weight', '.bin' 	! 电子权重
            call write_component_binary(filename, elec_local, nElec_local, 'weight', open_status, rank, MPI_COMM_WORLD, ierr)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
	if (rank == 0) print *, 'Done'
end subroutine write_particle_data
!*********************************************************************
subroutine write_component_binary(filename, p_array, n_local, component_name, open_status, rank, comm, ierr)	! 写入粒子数组分量
use pic_types, only: particle; use mpi
implicit none
    character(len=*), intent(in) :: filename
    type(particle), intent(in) :: p_array(:)	! 粒子数组指针
    integer, intent(in) :: n_local				! 本进程粒子数
    character(len=*), intent(in) :: component_name	! 要写入的分量
	character(len=*), intent(in) :: open_status
	integer, intent(in) :: rank, comm
	integer, intent(inout) :: ierr
    integer, parameter :: file_unit = 30
    integer :: i, io_stat, alloc_stat
	real(kind=8), allocatable :: data_buffer(:)
	
	if (n_local == 0) return
	allocate(data_buffer(n_local), stat=alloc_stat)		! 分配缓冲区
	if (alloc_stat /= 0) then
		print '("Rank ", I0, ": ERROR allocating data_buffer in write_component_binary. ", &
				"n_local = ", I0, ", component = ", A)', rank, n_local, trim(component_name)
		call mpi_abort(comm, -1, ierr)
	end if

	! 复制数据到连续缓冲区
	select case (trim(component_name))
        case ('x'); do i = 1, n_local; data_buffer(i) = p_array(i)%x; end do
        case ('y'); do i = 1, n_local; data_buffer(i) = p_array(i)%y; end do
        case ('z'); do i = 1, n_local; data_buffer(i) = p_array(i)%z; end do
        case ('vx'); do i = 1, n_local; data_buffer(i) = p_array(i)%vx; end do
        case ('vy'); do i = 1, n_local; data_buffer(i) = p_array(i)%vy; end do
        case ('vz'); do i = 1, n_local; data_buffer(i) = p_array(i)%vz; end do
        case ('weight'); do i = 1, n_local; data_buffer(i) = p_array(i)%weight; end do
        case default
            print '("Rank",I0,":ERROR in write_component_binary. Unknown component:",A)',rank,trim(component_name)
    end select
	
	io_stat = 0
	if (trim(open_status) == 'append') then	
		open(unit=file_unit, file=trim(filename), status='unknown', form='unformatted', &
			access='stream', position='append', action='write', iostat=io_stat)
	else
		open(unit=file_unit, file=trim(filename), status='replace', form='unformatted', &
			access='stream', action='write', iostat=io_stat)
	end if
	
	write(file_unit, iostat=io_stat) data_buffer	! 块写入
	if (io_stat /= 0) then
		print '("Rank",I0,":FATAL I/O ERROR *writing* to file:",A,".IOSTAT=",I0)',rank,trim(filename),io_stat
		call mpi_abort(comm, io_stat, ierr)
	end if
	close(file_unit)
	deallocate(data_buffer)
end subroutine write_component_binary
!*********************************************************************
subroutine dump_dynamic_binary(part, nLocal, species_in, step)	! 输出粒子位置(AoS 布局)
use mpi; use mesh_parameter; use pic_types
implicit none
	type(particle), intent(in) :: part(:)	! 本进程粒子数组
	integer, intent(in) :: nLocal
	character(len=*), intent(in) :: species_in
	integer, intent(in) :: step
	integer :: ierr, rank, nprocs, p
	character(len=64) :: fname_pos
	character(len=12) :: species
	real(kind=8), allocatable :: data_to_write(:,:)
	integer :: i
	integer :: active_p_count
	
	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
	
	species = adjustl(species_in)
	write(fname_pos, '("init_particles/",A,"_init_positon_",I0,".bin")') trim(species), step
	
	! 轮流写入
	do p=0, nprocs-1
		call MPI_Barrier(MPI_COMM_WORLD, ierr)
		if (rank == p) then
			active_p_count = 0
			do i=1, nLocal
				if(part(i)%weight /= 0.0d0) then		! 粒子权重不为0时输出
					active_p_count = active_p_count + 1
				end if
			end do
			
			if (active_p_count > 0) then
				allocate(data_to_write(3, active_p_count))
				active_p_count = 0
				do i=1, nLocal
					if (part(i)%weight /= 0.0d0) then		! 用于图像化，只写入权重非0的粒子
						active_p_count = active_p_count + 1
						data_to_write(1, active_p_count) = part(i)%x
						data_to_write(2, active_p_count) = part(i)%y
                        data_to_write(3, active_p_count) = part(i)%z
					end if
				end do
				
				open(unit=20, file=fname_pos, status='Unknown', action='write', access='stream', position='append')
					write(20) data_to_write
				close(20)
				deallocate(data_to_write)
			end if
		end if
	end do
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (rank == 0) print *, 'Done'
end subroutine dump_dynamic_binary
!*********************************************************************
end program pic_equilibrium
