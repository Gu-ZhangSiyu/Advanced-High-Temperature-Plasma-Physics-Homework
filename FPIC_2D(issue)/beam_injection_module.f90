!*********************************************************************
!   离子束流注入模块
!
! 空间分布	沿z轴方向注入，注入平面z坐标为z_min+21*dz
!			横截面采样： 在xy平面上采用均匀圆盘采样，离子在束斑圆内分布均匀
!
! 角度扰动	基于Box-Muller变换的二维独立高斯分布，束流的宏观几何发射度
! 速度分布	初始速度基于20keV能量 | 热运动特征速度基于0.5keV计算
!			热散射模型: 采用Box-Muller变换生成高斯分布的热速度
!			vx分量: 主速度xz投影 + 热速度扰动
!			vy分量:		0		 热速度扰动
!			vz分量: 主速度xz投影 +	 热速度扰动
!
! 调用说明	本模块支持入射中心调整
!			本模块兼容全同步与子循环异步两种运行模式
!			1. 全同步退化运行 (sub_ratio = 1)：
!			   将 dt_inject 与 dt_sim_base 均传递为主循环基准步长 dt。
!			2. 子循环异步运行 (sub_ratio > 1)：
!			   将 dt_inject 传递为离子大步长 dt_ion (即 dt * sub_ratio)；
!			   将 dt_sim_base 传递为主循环小步长 dt。
!			   此时，必须将该函数的调用移入离子触发分支(如 mod(idyn, sub_ratio) == 1)内。
!*********************************************************************
module beam_injection_module
    use mesh_parameter, only: dz, z_min, y_plane
    use physicals_parameter, only: q, M_ION, pi
    use pic_types, only: particle
	use control_parameter
    use Domain_Decomposition_Module, only: domain_decomposition, Owner_Of_Position
    use Particle_Migration_Module, only: Ensure_Particle_Capacity
    implicit none
contains
!*********************************************************************
subroutine Inject_Beam_Species(part_array, nP_active, nP_max, dt_inject, dt_sim_base, idyn, domain, &
                               m_species, q_species, total_target_particles, macro_weight, &
                               E_beam_keV, global_injected_total, fractional_particles) 
                               
    type(particle), pointer, intent(inout) :: part_array(:)
    integer :: inject_count, i, new_idx, source_rank
    integer, intent(in) :: idyn
    integer, intent(inout) :: nP_max
    integer, intent(inout) :: nP_active
    type(domain_decomposition), intent(in) :: domain
    integer(kind=8), intent(inout) :: global_injected_total
    real(8), intent(in) :: dt_inject, dt_sim_base, m_species, q_species, total_target_particles, macro_weight, E_beam_keV
    real(8), intent(inout) :: fractional_particles
    real(8) :: current_time, macro_per_step, v0, v_thermal, theta_pitch, rand_r, rand_u1, rand_v1, rand_u2, rand_v2, &
               E_joule, delta_theta, delta_phi, v0_x, v0_y, v0_z
    real(8), dimension(3) :: v_vec, pos_vec

    current_time = dble(idyn) * dt_sim_base
    if (current_time > injection_duration) return

    ! 注入平面只属于一个空间子区域，只有该 rank 生成粒子。
    source_rank = Owner_Of_Position(z_min + inject_z_offset * dz, domain)
    if (domain%rank /= source_rank) return
    macro_per_step = total_target_particles * dt_inject / injection_duration
    fractional_particles = fractional_particles + macro_per_step
    inject_count = floor(fractional_particles)
    fractional_particles = fractional_particles - dble(inject_count)
    if (inject_count == 0) return

    if (inject_count > 0) then
        call Ensure_Particle_Capacity(part_array, nP_active, nP_max, nP_active + inject_count)
    end if

    ! 基础速度与热散布计算
    E_joule = E_beam_keV * 1000.0d0 * abs(q_species)
    v0 = sqrt(2.0d0 * E_joule / m_species)
    v_thermal = sqrt(2.0d0 * (E_thermal_eV * abs(q_species)) / m_species)
    theta_pitch = pitch_angle_deg * pi / 180.0d0

    do i = 1, inject_count
        new_idx = nP_active + 1

        ! 2D XZ 空间位置采样：圆盘横截面退化为 x 方向束斑线段。
        call random_number(rand_u1)
		rand_r = beam_spot_radius * (2.0d0 * rand_u1 - 1.0d0)
		pos_vec(1) = beam_center_x + rand_r
		pos_vec(2) = y_plane
		pos_vec(3) = z_min + inject_z_offset * dz

        ! 速度采样准备
        call random_number(rand_u1); call random_number(rand_v1); call random_number(rand_u2); call random_number(rand_v2)
        if (rand_u1 < 1.0d-12) rand_u1 = 1.0d-12; if (rand_u2 < 1.0d-12) rand_u2 = 1.0d-12

        ! 施加几何发散
        delta_theta = (angle_perturb_deg * pi / 180.0d0) * sqrt(-2.0d0*log(rand_u1)) * cos(2.0d0*pi*rand_v1)
        delta_phi = (angle_perturb_deg * pi / 180.0d0) * sqrt(-2.0d0*log(rand_u1)) * sin(2.0d0*pi*rand_v1)
        v0_x = v0 * sin(theta_pitch + delta_theta); v0_y = v0 * sin(delta_phi); v0_z = v0 * cos(theta_pitch + delta_theta)

        ! 叠加三维热扰动
        v_vec(1) = v0_x + v_thermal * sqrt(-2.0d0*log(rand_u2)) * cos(2.0d0*pi*rand_v2)
        v_vec(2) = v0_y + v_thermal * sqrt(-2.0d0*log(rand_u2)) * sin(2.0d0*pi*rand_v2)

        call random_number(rand_u1); call random_number(rand_v1)
        if (rand_u1 < 1.0d-12) rand_u1 = 1.0d-12
        v_vec(3) = v0_z + v_thermal * sqrt(-2.0d0*log(rand_u1)) * cos(2.0d0*pi*rand_v1)
        
        ! 粒子赋值
        part_array(new_idx)%x = pos_vec(1); part_array(new_idx)%y = pos_vec(2); part_array(new_idx)%z = pos_vec(3)
        part_array(new_idx)%vx = v_vec(1); part_array(new_idx)%vy = v_vec(2); part_array(new_idx)%vz = v_vec(3)
        part_array(new_idx)%weight = macro_weight
        part_array(new_idx)%id = global_injected_total + int(i, 8)
        nP_active = new_idx
    end do
    global_injected_total = global_injected_total + int(inject_count, 8)
end subroutine Inject_Beam_Species
!*********************************************************************
end module beam_injection_module
