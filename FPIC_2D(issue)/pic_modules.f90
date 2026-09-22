!*********************************************************************
! 参数设置
!*********************************************************************
module mesh_parameter
    implicit none
    integer, parameter :: nx = 129, ny = 1, nz = 129		! 网格数 (x,y,z)
    real(8), parameter :: x_min = -2.0d-2		! 网格的物理边界 (m)
    real(8), parameter :: x_max =  2.0d-2
    real(8), parameter :: y_min = -2.0d-2
    real(8), parameter :: y_max =  2.0d-2
    real(8), parameter :: y_plane = 0.0d0
    real(8), parameter :: z_min = -2.0d-2
    real(8), parameter :: z_max =  2.0d-2
    real(8), parameter :: dx = (x_max - x_min) / dble(nx - 1)	! 网格间距(均匀网格)
    real(8), parameter :: dy = y_max - y_min
    real(8), parameter :: dz = (z_max - z_min) / dble(nz - 1)	
	real(8), parameter :: inv_dx = 1.0d0 / dx
	real(8), parameter :: inv_dy = 1.0d0 / dy
	real(8), parameter :: inv_dz = 1.0d0 / dz
end module mesh_parameter
!*********************************************************************
module physicals_parameter
    implicit none
	real(8), parameter :: pi = 3.141592653589793d0	! 圆周率
	real(8), parameter :: q = 1.60217733d-19		! 元电荷 (C)
	real(8), parameter :: c = 2.99792458d8			! 光速
	real(8), parameter :: inv_eps0 = 1.0d0 / 8.8541878128d-12
	real(8), parameter :: M_ELEC = 9.1093837139d-31	! 电子质量
	real(8), parameter :: M_ION = 25.0d0 * M_ELEC	! 离子质量（缩放质量比）
	real(8), parameter :: c2 = c*c
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
	
	! 负载均衡
    logical, parameter :: enable_load_balance = .false.			! 动态负载均衡开关
    integer, parameter :: load_balance_interval = 1000 			! 两次运行中检查至少间隔的电子步数
    real(8), parameter :: load_balance_threshold = 1.15d0 		! 旧分区最大预计成本
    real(8), parameter :: load_balance_min_gain = 0.05d0 		! 新分区最大预计成本至少降低5%才实施
    real(8), parameter :: load_balance_cell_cost = 0.1d0		! 每网格节点相对成本
    real(8), parameter :: load_balance_electron_cost = 1.0d0	! 每电子宏粒子相对成本
    real(8), parameter :: load_balance_ion_cost = 1.0d0 		! 每离子宏粒子成本，除以sub_ratio
    real(8), parameter :: load_balance_max_slab_factor = 2.0d0 	! 单rank最大网格层数相对初始均匀层数上限的倍数

	! 功能开关
    integer, parameter :: sub_ratio = 10						! 电子离子异步倍率（设置为1同步）
    logical, parameter :: load_background_particles = .true. 	! 初始读取背景粒子（重启时忽略）
    logical, parameter :: enable_beam_injection = .false.		! 束流注入开关
    logical, parameter :: enable_RF = .false.					! 外加RF电场开关
	
    ! 粒子权重设置
	! 注: 手动设置的背景粒子权重与入射粒子权重无关
    real(8), parameter :: background_weight_ion = 1.0d0  		! 手动设置背景离子权重
    real(8), parameter :: background_weight_elec = 1.0d0 		! 手动设置背景电子权重
    logical, parameter :: background_weights_from_file = .true.	! 2.pic_init.f90 输出逐粒子权重
    real(8), parameter :: legacy_restart_weight_ion = 1.25d3		! 无动态权重文件时重启使用的离子统一权重
    real(8), parameter :: legacy_restart_weight_elec = 1.25d3		! 无动态权重文件时重启使用的电子统一权重
    integer(kind=8), parameter :: nIon  = 10000000					! 初始超离子数
    integer(kind=8), parameter :: nElec = 10000000					! 初始超电子数
    integer :: nStepsDyn = 1000								! 主循环步数
    real(8), parameter :: dt = 5d-13							! 时间步长(s)
	real(8), parameter :: inv_dt = 1.0d0 / dt
    integer :: restart_step = 0									! 0:初始运行，非0:从该步读取文件重启

    ! 输出控制
    integer :: freq_ion = 10				! 离子输出频率
    integer :: freq_elec = 10				! 电子输出频率
    integer :: freq_field = 10			! 输出场频率
    integer :: sample_rate_ion = 1			! 离子降采样保存率(每隔多少个ID记录一个粒子)
    integer :: sample_rate_elec = 1			! 电子降采样保存率
    logical :: out_ion = .true.				! 是否输出离子(只影响文件保存)
    logical :: out_elec = .true.			! 是否输出电子
    logical :: out_field = .true.			! 是否输出场
	
    ! 束流注入控制参数
    real(8), parameter :: E_beam_ion_keV = 10.0d0       ! 离子束能量 (keV)
    real(8), parameter :: E_beam_elec_keV = 5.4461d-3    ! 电子束能量 (keV)
    real(8), parameter :: total_target_ions = 2.0d7  	! 超离子总数
    real(8), parameter :: total_target_elec = 2.0d7  	! 超电子总数
    real(8), parameter :: weight_ion = 1.25d3            ! 超离子权重
    real(8), parameter :: weight_elec = 1.25d3           ! 超电子权重
	real(8), parameter :: pitch_angle_deg = 0.0d0      ! 螺旋入射角 (度)
    real(8), parameter :: beam_spot_radius = 0.005d0    ! 束斑半径 (m)
    real(8), parameter :: inject_z_offset = 0.0d0      ! z轴入射平面偏移量 (网格)
    real(8), parameter :: injection_duration = 4.0d-7   ! 注入持续时间 (s)
    real(8), parameter :: E_thermal_eV = 10.0d0         ! 热速度对应能量 (eV)
    real(8), parameter :: angle_perturb_deg = 0.5d0     ! 几何角度扰动标准差 (度)
	real(8), parameter :: beam_center_x = 0.0d0  		! 2D入射位置X方向偏移 (m)

    ! 外加RF电场参数
    real(8), parameter :: CARA_RF_phase = 0.0d0              ! RF初始相位 (rad)
    real(8), parameter :: CARA_RF_Ey_amplitude = 2.0d5       ! RF电场峰值 (V/m)
    real(8), parameter :: CARA_RF_omega = 2.394709329d8       ! RF角频率 (rad/s), 对应38.1129827 MHz
end module control_parameter
!*********************************************************************
