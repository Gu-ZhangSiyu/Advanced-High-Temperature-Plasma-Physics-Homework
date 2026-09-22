program generate_3d_reconnection_equilibrium
    implicit none
    integer :: i, j, k, ic, jc, kc
    integer, parameter :: nx = 257, ny = 257, nz = 3
    real(8), parameter :: dx = 0.5d-3, dy = 0.5d-3, dz = 1.0d-3     ! 单网格尺寸
    real(8), parameter :: e_charge = 1.60217663d-19
    real(8), parameter :: k_B = 1.380649d-23
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: Te_eV = 1000.0d0
    real(8), parameter :: Ti_eV = 1000.0d0
    real(8), parameter :: n0 = 1.0d17          ! 中心片峰值密度
    real(8), parameter :: n_bg = 2.0d16        ! 维持计算稳定的本底密度
    real(8), parameter :: B0 = 0.05d0          ! 重联反转磁场
    real(8), parameter :: lambda = 24.0d0 * dy ! 电流片半宽度 (约0.5个电子趋肤深度，激进一点以加速重联)
    real(8), parameter :: pert = 0.05d0        ! 初始磁扰动振幅 (5% B0)
    real(8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
    real(8), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
    real(8), allocatable :: dens(:,:,:), prs(:,:,:), psi(:,:,:)
    real(8) :: x, y, y_c, x_c, Total_T, d_psi
    character(len=20) :: fmt_xy, fmt_xz, fmt_yz
    
    allocate(Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz))
    allocate(Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz))
    allocate(dens(nx,ny,nz), prs(nx,ny,nz), psi(nx,ny,nz))
    
    Total_T = (Te_eV + Ti_eV) * e_charge
    y_c = dble(ny)/2.0d0 * dy  ! Y方向中心面 (重联发生在此处)
    x_c = dble(nx)/2.0d0 * dx  ! X方向中心面 (X点初始位置)
    ic = nx / 2
    jc = ny / 2
    kc = nz / 2
    
    do k = 1, nz
        do j = 1, ny
            y = (dble(j) - 0.5d0) * dy
            
            do i = 1, nx
                x = (dble(i) - 0.5d0) * dx
                
                ! 基础哈里斯片磁场 (Bx 反转) 及 引导磁场 (Bz)
                Bx(i,j,k) = B0 * tanh((y - y_c) / lambda)
                By(i,j,k) = 0.0d0
                Bz(i,j,k) = 0.0d0
                
                ! 基础磁通量函数 psi (这里实质上是 Z 方向的矢量势 A_z)
                psi(i,j,k) = -B0 * lambda * log(cosh((y - y_c) / lambda))
                
                ! X-Y 平面施加重联微扰 (X点位于域中心)
                d_psi = -pert * B0 * lambda * cos(2.0d0 * pi * (x - x_c) / (dble(nx)*dx)) * exp(-((y - y_c) / lambda)**2)
                psi(i,j,k) = psi(i,j,k) + d_psi
                
                ! 扰动磁场计算 (同步将对应的符号翻转)
                if (j > 1 .and. j < ny .and. i > 1 .and. i < nx) then
                    Bx(i,j,k) = Bx(i,j,k) + (pert * B0 * lambda * cos(2.0d0 * pi * (x - x_c) / (dble(nx)*dx)) * &
                                (-2.0d0 * (y - y_c) / lambda**2) * exp(-((y - y_c) / lambda)**2))
                    By(i,j,k) = pert * B0 * lambda * (2.0d0 * pi / (dble(nx)*dx)) * &
                                sin(2.0d0 * pi * (x - x_c) / (dble(nx)*dx)) * exp(-((y - y_c) / lambda)**2)
                end if
                
                ! 密度剖面与压力闭合
                dens(i,j,k) = n0 * (1.0d0 / cosh((y - y_c) / lambda))**2 + n_bg
                prs(i,j,k) = dens(i,j,k) * Total_T
                
                ! 初始电场
                Ex(i,j,k) = 0.0d0; Ey(i,j,k) = 0.0d0; Ez(i,j,k) = 0.0d0
            end do
        end do
    end do
    
    open(unit=10, file='out_3d/p_3d.bin', form='unformatted', access='stream', status='replace'); write(10) prs; close(10)
    open(unit=10, file='out_3d/psi_3d.bin', form='unformatted', access='stream', status='replace'); write(10) psi; close(10)
    open(unit=10, file='out_3d/den_3d.bin', form='unformatted', access='stream', status='replace'); write(10) dens; close(10)
    open(unit=10, file='out_3d/Bx_3d.bin', form='unformatted', access='stream', status='replace'); write(10) Bx; close(10)
    open(unit=10, file='out_3d/By_3d.bin', form='unformatted', access='stream', status='replace'); write(10) By; close(10)
    open(unit=10, file='out_3d/Bz_3d.bin', form='unformatted', access='stream', status='replace'); write(10) Bz; close(10)
    open(unit=10, file='out_3d/Ex_3d.bin', form='unformatted', access='stream', status='replace'); write(10) Ex; close(10)
    open(unit=10, file='out_3d/Ey_3d.bin', form='unformatted', access='stream', status='replace'); write(10) Ey; close(10)
    open(unit=10, file='out_3d/Ez_3d.bin', form='unformatted', access='stream', status='replace'); write(10) Ez; close(10)
    write(fmt_xy, '("(",I0,"E15.6)")') nx
    write(fmt_xz, '("(",I0,"E15.6)")') nx
    write(fmt_yz, '("(",I0,"E15.6)")') ny
    open(unit=50, file='out_2d/bx_xy.txt', status='replace')
    open(unit=51, file='out_2d/by_xy.txt', status='replace')
    open(unit=52, file='out_2d/bz_xy.txt', status='replace')
    open(unit=53, file='out_2d/dens_xy.txt', status='replace')
    open(unit=54, file='out_2d/psi_xy.txt', status='replace')
    open(unit=55, file='out_2d/prs_xy.txt', status='replace')
    do j = 1, ny
        write(50, fmt_xy) Bx(:, j, kc)
        write(51, fmt_xy) By(:, j, kc)
        write(52, fmt_xy) Bz(:, j, kc)
        write(53, fmt_xy) dens(:, j, kc)
        write(54, fmt_xy) psi(:, j, kc)
        write(55, fmt_xy) prs(:, j, kc)
    end do
    close(50); close(51); close(52); close(53); close(54); close(55)
end program generate_3d_reconnection_equilibrium