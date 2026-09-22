!***************************************************************
!	Coded By Siyu Zhang     November. 30, 2025
!
!   Boris粒子推进数学求解器（兼容相对论形式）适配Yee网格版
!
!	边界处理:
!	boundary_control=1: 报告模式: 开启飞出模拟区域检测，如果飞出边界会给status设置一个非0值后直接return
!	boundary_control=0: 穿越模式: 开启穿越模式的边界判断，已经飞出边界的粒子失去场的束缚
!	
!	相对论粒子：
!	relativistic_control=1: 相对论粒子
!	relativistic_control=0: 非相对论粒子
!
!	调用示例:
!	call Push_Particle(pos_in=my_pos_in, vel_in=my_vel_in, &
!						qm=qm0, dt=dt, Ex=Ex, Ey=Ey, Ez=Ez, Bx=Bx, By=By, Bz=Bz, &
!                       grid_origin=origin_coords, grid_steps=cell_size, grid_dims=[nx, ny, nz], &
!                       pos_out=my_pos_out, vel_out=my_vel_out, &
!                       status=particle_status, boundary_control=1, relativistic_control=1)
! 	
!	传递值说明:
!	pos_in			粒子初始位置						数组，例如 my_pos(3) = [1.0, 2.5, 3.0]
!	vel_in			粒子的初始速度						数组，例如 my_vel(3) = [1e5, -2e5, 0.0]
!	qm				粒子的荷质比(q/m)					标量值，例如 cm = -1.7588e11
!	qm（相对论）		电荷与粒子静止质量的比值(q/m0)		标量值，例如 cm0 = -1.7588e11
!	dt				时间步长 							标量值，例如 time_step = 1.0e-11		
!	Ex, Ey, Ez		电场分量							三个三维数组，例如 E_field_x(nx, ny, nz)
!	Bx, By, Bz		磁场分量							三个三维数组，例如 B_field_x(nx, ny, nz)
!	grid_origin		网格原点(索引1,1,1处)物理坐标		大小为3的数组，例如 origin_coords = [0.0, 0.0, 0.0]
!	grid_steps		网格在三个方向的步长(dx,dy,dz)		大小为3的数组，例如 cell_size = [0.01, 0.01, 0.02]	
!	grid_dims		网格在三个方向的维度(nx,ny,nz)		大小为3的数组，例如用[nx, ny, nz]构造
!	pos_out			经过 Delta_t 后粒子的新位置		提供一个已声明的大小为3的数组，用于接收结果
!	vel_out			经过 Delta_t 后粒子的新速度		提供一个已声明的大小为3的数组，用于接收结果
!	status			粒子是否越界的状态码				提供一个已声明的integer变量，用于接收状态码
!	boundary_control		边界处理模式开关			布尔值
!	relativistic_control	粒子相对论模式开关			布尔值
!***************************************************************
module particle_solver_yee
    implicit none
    private                     ! 默认模块内所有内容为私有
    public :: Push_Particle     ! 只将核心求解器公开

contains

    subroutine Push_Particle(pos_in, vel_in, qm, dt, Ex, Ey, Ez, Bx, By, Bz, &
                            grid_origin, grid_steps, grid_dims, pos_out, vel_out, &
                            status, boundary_control, relativistic_control)
        implicit none

        ! input
        real(8), intent(in) :: pos_in(3)      ! 粒子初始位置 [x, y, z]
        real(8), intent(in) :: vel_in(3)      ! 粒子初始速度 [vx, vy, vz]
        real(8), intent(in) :: qm             ! 荷质比 (q / m0)
        real(8), intent(in) :: dt             ! 时间步长
        real(8), intent(in) :: Ex(0:,0:,0:), Ey(0:,0:,0:), Ez(0:,0:,0:)
        real(8), intent(in) :: Bx(0:,0:,0:), By(0:,0:,0:), Bz(0:,0:,0:)
        real(8), intent(in) :: grid_origin(3) ! 网格原点(1,1,1)对应的物理坐标 [x0, y0, z0]
        real(8), intent(in) :: grid_steps(3)  ! 网格步长 [dx, dy, dz]
        integer, intent(in) :: grid_dims(3)   ! 网格维度 [nx, ny, nz]
        integer, intent(in), optional :: boundary_control ! 控制边界检查: 0=钳位(默认), 1=报告
        integer, intent(in), optional :: relativistic_control ! 控制计算模式: 0=非相对论(默认), 1=相对论

        ! output
        real(8), intent(out) :: pos_out(3)    ! 粒子在 dt 之后的新位置 [x, y, z]
        real(8), intent(out) :: vel_out(3)    ! 粒子在 dt 之后的新速度 [vx, vy, vz]
        integer, intent(out) :: status        ! 粒子状态: 0=界内, +/-1=X方向越界, +/-2=Y方向越界, +/-3=Z方向越界

        ! for cic
        integer :: mode, rel_mode
        real(8) :: logic_pos(3)
        integer :: i_n, j_n, k_n      ! 整数/对齐网格索引
        integer :: i_h, j_h, k_h      ! 半整数/偏移网格索引
        integer :: i0, j0, k0
        real(8) :: wxn_0, wxn_1, wyn_0, wyn_1, wzn_0, wzn_1 ! 整数网格权重
        real(8) :: wxh_0, wxh_1, wyh_0, wyh_1, wzh_0, wzh_1 ! 半整数网格权重
        real(8) :: shift_pos
        real(8) :: E_loc(3), B_loc(3)

		! for boris (non-relativistic)
        real(8) :: v_minus(3), v_prime(3), v_plus(3)
        real(8) :: t_sq

        ! for boris (relativistic)
        real(8), parameter :: C_LIGHT = 299792458.0d0
        real(8), parameter :: C_LIGHT_SQ = C_LIGHT**2
        real(8) :: u_in(3), u_minus(3), u_prime(3), u_plus(3), u_out(3)
        real(8) :: gamma
        real(8) :: vel_sq

        ! for boris (common)
        real(8) :: half_dt
        real(8) :: t(3), s(3)

        ! 确定边界处理模式
        if (present(boundary_control)) then
            mode = boundary_control
        else
            mode = 0
        end if

        ! 确定计算模式
        if (present(relativistic_control)) then
            rel_mode = relativistic_control
        else
            rel_mode = 0
        end if

        ! 计算粒子在网格中的逻辑坐标 (从0开始计数)
        logic_pos(1) = (pos_in(1) - grid_origin(1)) / grid_steps(1)
        logic_pos(2) = (pos_in(2) - grid_origin(2)) / grid_steps(2)
        logic_pos(3) = (pos_in(3) - grid_origin(3)) / grid_steps(3)
        
        ! 找到包围粒子的立方体单元的左下角索引
        i0 = int(logic_pos(1)) + 1
        j0 = int(logic_pos(2)) + 1
        k0 = int(logic_pos(3)) + 1

        ! 边界检查
        status = 0
        if (mode == 1) then
            if (i0 < 1) status = -1
            if (i0 >= grid_dims(1)) status = 1
            if (j0 < 1) status = -2
            if (j0 >= grid_dims(2)) status = 2
            if (k0 < 1) status = -3
            if (k0 >= grid_dims(3)) status = 3

            ! 粒子出界直接返回
            if (status /= 0) then
                pos_out = pos_in
                vel_out = vel_in
                return
            end if
        end if

        ! 内存安全网
        if (mode == 0) then
            if (i0 < 0 .or. i0 > grid_dims(1) + 1 .or. &
                j0 < 0 .or. j0 > grid_dims(2) + 1 .or. &
                k0 < 0 .or. k0 > grid_dims(3) + 1) then
                status = 9    ! 严重越界码
                pos_out = pos_in; vel_out = vel_in
                return
            end if
        end if

        ! Yee Grid 插值
        ! 计算整数/对齐网格的索引和权重
        i_n = floor(logic_pos(1)) + 1
        wxn_1 = logic_pos(1) - real(i_n - 1, 8)
        wxn_0 = 1.0d0 - wxn_1
        
        j_n = floor(logic_pos(2)) + 1
        wyn_1 = logic_pos(2) - real(j_n - 1, 8)
        wyn_0 = 1.0d0 - wyn_1
        
        k_n = floor(logic_pos(3)) + 1
        wzn_1 = logic_pos(3) - real(k_n - 1, 8)
        wzn_0 = 1.0d0 - wzn_1

        ! 计算半整数/偏移网格的索引和权重
        shift_pos = logic_pos(1) - 0.5d0
        i_h = floor(shift_pos) + 1
        wxh_1 = shift_pos - real(i_h - 1, 8)
        wxh_0 = 1.0d0 - wxh_1

        shift_pos = logic_pos(2) - 0.5d0
        j_h = floor(shift_pos) + 1
        wyh_1 = shift_pos - real(j_h - 1, 8)
        wyh_0 = 1.0d0 - wyh_1

        shift_pos = logic_pos(3) - 0.5d0
        k_h = floor(shift_pos) + 1
        wzh_1 = shift_pos - real(k_h - 1, 8)
        wzh_0 = 1.0d0 - wzh_1

        ! 组装场分量
        E_loc(1) = (Ex(i_h, j_n, k_n)*wyn_0*wzn_0 + Ex(i_h, j_n+1, k_n)*wyn_1*wzn_0 + &
                    Ex(i_h, j_n, k_n+1)*wyn_0*wzn_1 + Ex(i_h, j_n+1, k_n+1)*wyn_1*wzn_1) * wxh_0 + &
                    (Ex(i_h+1, j_n, k_n)*wyn_0*wzn_0 + Ex(i_h+1, j_n+1, k_n)*wyn_1*wzn_0 + &
                    Ex(i_h+1, j_n, k_n+1)*wyn_0*wzn_1 + Ex(i_h+1, j_n+1, k_n+1)*wyn_1*wzn_1) * wxh_1

        E_loc(2) = (Ey(i_n, j_h, k_n)*wxn_0*wzn_0 + Ey(i_n+1, j_h, k_n)*wxn_1*wzn_0 + &
                    Ey(i_n, j_h, k_n+1)*wxn_0*wzn_1 + Ey(i_n+1, j_h, k_n+1)*wxn_1*wzn_1) * wyh_0 + &
                    (Ey(i_n, j_h+1, k_n)*wxn_0*wzn_0 + Ey(i_n+1, j_h+1, k_n)*wxn_1*wzn_0 + &
                    Ey(i_n, j_h+1, k_n+1)*wxn_0*wzn_1 + Ey(i_n+1, j_h+1, k_n+1)*wxn_1*wzn_1) * wyh_1

        E_loc(3) = (Ez(i_n, j_n, k_h)*wxn_0*wyn_0 + Ez(i_n+1, j_n, k_h)*wxn_1*wyn_0 + &
                    Ez(i_n, j_n+1, k_h)*wxn_0*wyn_1 + Ez(i_n+1, j_n+1, k_h)*wxn_1*wyn_1) * wzh_0 + &
                    (Ez(i_n, j_n, k_h+1)*wxn_0*wyn_0 + Ez(i_n+1, j_n, k_h+1)*wxn_1*wyn_0 + &
                    Ez(i_n, j_n+1, k_h+1)*wxn_0*wyn_1 + Ez(i_n+1, j_n+1, k_h+1)*wxn_1*wyn_1) * wzh_1

        B_loc(1) = (Bx(i_n, j_h, k_h)*wyh_0*wzh_0 + Bx(i_n, j_h+1, k_h)*wyh_1*wzh_0 + &
                    Bx(i_n, j_h, k_h+1)*wyh_0*wzh_1 + Bx(i_n, j_h+1, k_h+1)*wyh_1*wzh_1) * wxn_0 + &
                    (Bx(i_n+1, j_h, k_h)*wyh_0*wzh_0 + Bx(i_n+1, j_h+1, k_h)*wyh_1*wzh_0 + &
                    Bx(i_n+1, j_h, k_h+1)*wyh_0*wzh_1 + Bx(i_n+1, j_h+1, k_h+1)*wyh_1*wzh_1) * wxn_1

        B_loc(2) = (By(i_h, j_n, k_h)*wxh_0*wzh_0 + By(i_h,   j_n, k_h+1)*wxh_0*wzh_1 + &
                    By(i_h+1, j_n, k_h)*wxh_1*wzh_0 + By(i_h+1, j_n, k_h+1)*wxh_1*wzh_1) * wyn_0 + &
                    (By(i_h, j_n+1, k_h)*wxh_0*wzh_0 + By(i_h,   j_n+1, k_h+1)*wxh_0*wzh_1 + &
                    By(i_h+1, j_n+1, k_h)*wxh_1*wzh_0 + By(i_h+1, j_n+1, k_h+1)*wxh_1*wzh_1) * wyn_1

        B_loc(3) = (Bz(i_h, j_h, k_n)*wxh_0*wyh_0 + Bz(i_h+1, j_h, k_n)*wxh_1*wyh_0 + &
                    Bz(i_h, j_h+1, k_n)*wxh_0*wyh_1 + Bz(i_h+1, j_h+1, k_n)*wxh_1*wyh_1) * wzn_0 + &
                    (Bz(i_h, j_h, k_n+1)*wxh_0*wyh_0 + Bz(i_h+1, j_h, k_n+1)*wxh_1*wyh_0 + &
                    Bz(i_h, j_h+1, k_n+1)*wxh_0*wyh_1 + Bz(i_h+1, j_h+1, k_n+1)*wxh_1*wyh_1) * wzn_1

        ! Boris推进
        half_dt = 0.5d0 * dt
        if (rel_mode == 1) then
            ! Relativistic Boris推进
            vel_sq = vel_in(1)**2 + vel_in(2)**2 + vel_in(3)**2		! 初速度转换为相对论性动量
            if (vel_sq >= C_LIGHT_SQ) then
                vel_sq = C_LIGHT_SQ * 0.999999999999999d0	! 非物理速度限制在c以下
            end if
            gamma = 1.0d0 / sqrt(1.0d0 - vel_sq / C_LIGHT_SQ)
            u_in = gamma * vel_in
            u_minus = u_in + qm * E_loc * half_dt		! 电场作用下的前半步推进

            ! 磁场作用下的旋转
            gamma = sqrt(1.0d0 + (u_minus(1)**2 + u_minus(2)**2 + u_minus(3)**2) / C_LIGHT_SQ)
            t = qm * B_loc * half_dt / gamma
            s = 2.0d0 * t / (1.0d0 + (t(1)**2 + t(2)**2 + t(3)**2))
            u_prime = u_minus + cross_product(u_minus, t)
            u_plus = u_minus + cross_product(u_prime, s)

            u_out = u_plus + qm * E_loc * half_dt		! 电场作用下的后半步推进

            ! 将相对论性动量转换回速度
            gamma = sqrt(1.0d0 + (u_out(1)**2 + u_out(2)**2 + u_out(3)**2) / C_LIGHT_SQ)
            vel_out = u_out / gamma
        else
            ! Non-relativistic Boris推进
            v_minus = vel_in + qm * E_loc * half_dt			! 电场作用下的前半步推进
            t = qm * B_loc * half_dt						! 磁场作用下的旋转
            t_sq = t(1)**2 + t(2)**2 + t(3)**2
            s = 2.0d0 * t / (1.0d0 + t_sq)					! 计算辅助旋转矢量，避免直接算三角函数
            
            ! 通过两次叉乘完成旋转
            v_prime = v_minus + cross_product(v_minus, t)	
            v_plus  = v_minus + cross_product(v_prime, s)
            
            vel_out = v_plus + qm * E_loc * half_dt			! 电场作用下的后半步推进
        end if

        pos_out = pos_in + vel_out * dt					! 更新粒子位置
        
    end subroutine Push_Particle

    ! 计算三维矢量叉乘
    pure function cross_product(a, b) result(c)
        real(8), intent(in) :: a(3), b(3)
        real(8) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

end module particle_solver_yee