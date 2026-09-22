! 2D XZ Harris 重联背景场生成器
!
! 原始测试程序使用 (nx,ny,nz)=(257,257,3)，物理平面为 XY，z=3 是薄方向。
! 当前 Full-PIC 工程统一使用 XZ 平面，因此这里把原来的 y 方向旋转为 z，
! 并将薄方向压缩为 ny=1。输出路径保持与主程序 data_loader_module 一致：
!     out_xz/{Bx,By,Bz,Ex,Ey,Ez}_xz.txt
!
program generate_2d_xz_harris
    implicit none
    integer, parameter :: nx = 129, ny = 1, nz = 129
    real(8), parameter :: x_min = -2.0d-2, x_max = 2.0d-2
    real(8), parameter :: y_min = -2.0d-2, y_max = 2.0d-2
    real(8), parameter :: z_min = -2.0d-2, z_max = 2.0d-2
    real(8), parameter :: e_charge = 1.60217663d-19
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: Te_eV = 1000.0d0
    real(8), parameter :: Ti_eV = 1000.0d0
    real(8), parameter :: n0 = 1.0d17
    real(8), parameter :: n_bg = 2.0d16
    real(8), parameter :: B0 = 0.05d0
    real(8), parameter :: pert = 0.05d0
    real(8), parameter :: dx = (x_max - x_min) / dble(nx - 1)
    real(8), parameter :: dz = (z_max - z_min) / dble(nz - 1)
    real(8), parameter :: y_depth = y_max - y_min
    real(8), parameter :: lambda = 24.0d0 * dz
    real(8), parameter :: x_center = 0.0d0, z_center = 0.0d0

    real(8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:)
    real(8), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
    real(8), allocatable :: dens(:,:,:), prs(:,:,:), psi(:,:,:)
    integer :: i, k
    real(8) :: xx, zz, s, phase, dpsi, total_temperature

    allocate(Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz))
    allocate(Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz))
    allocate(dens(nx,ny,nz), prs(nx,ny,nz), psi(nx,ny,nz))

    total_temperature = (Te_eV + Ti_eV) * e_charge
    do k = 1, nz
        zz = z_min + dble(k - 1) * dz
        s = (zz - z_center) / lambda
        do i = 1, nx
            xx = x_min + dble(i - 1) * dx
            phase = 2.0d0 * pi * (xx - x_center) / (x_max - x_min)

            ! A_y(x,z) 的 XZ 旋转版本：Bx=-dA_y/dz, Bz=dA_y/dx。
            dpsi = -pert * B0 * lambda * cos(phase) * exp(-s*s)
            psi(i,1,k) = -B0 * lambda * log(cosh(s)) + dpsi
            Bx(i,1,k) = B0 * tanh(s) - &
                        2.0d0 * pert * B0 * s * cos(phase) * exp(-s*s)
            By(i,1,k) = 0.0d0
            Bz(i,1,k) = pert * B0 * lambda * (2.0d0*pi/(x_max-x_min)) * &
                        sin(phase) * exp(-s*s)

            dens(i,1,k) = n0 / (cosh(s)**2) + n_bg
            prs(i,1,k) = dens(i,1,k) * total_temperature
            Ex(i,1,k) = 0.0d0
            Ey(i,1,k) = 0.0d0
            Ez(i,1,k) = 0.0d0
        end do
    end do

    call write_xz('out_xz/Bx_xz.txt', Bx)
    call write_xz('out_xz/By_xz.txt', By)
    call write_xz('out_xz/Bz_xz.txt', Bz)
    call write_xz('out_xz/Ex_xz.txt', Ex)
    call write_xz('out_xz/Ey_xz.txt', Ey)
    call write_xz('out_xz/Ez_xz.txt', Ez)
    call write_xz('out_xz/density_xz.txt', dens)
    call write_xz('out_xz/pressure_xz.txt', prs)
    call write_xz('out_xz/psi_xz.txt', psi)

    print '(A,I0,A,I0,A)', 'Generated 2D XZ Harris field: ', nx, ' x ', nz, ' nodes.'
    print '(A,ES12.5,A,ES12.5)', 'x range = [', x_min, ', ', x_max, ']'
    print '(A,ES12.5,A,ES12.5)', 'z range = [', z_min, ', ', z_max, ']'
    print '(A,ES12.5)', 'effective y depth = ', y_depth

contains
    subroutine write_xz(filename, data)
        character(len=*), intent(in) :: filename
        real(8), intent(in) :: data(:,:,:)
        integer :: unit_no, i, k, ios

        open(newunit=unit_no, file=filename, form='formatted', status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) error stop 'cannot open XZ field text file'
        ! 每行对应一个 z 网格层，行内按 x=1...nx 写出；不写表头。
        do k = 1, size(data, 3)
            write(unit_no, '(*(ES24.16,1X))', iostat=ios) (data(i,1,k), i=1,size(data,1))
            if (ios /= 0) error stop 'cannot write XZ field text file'
        end do
        close(unit_no)
    end subroutine write_xz
end program generate_2d_xz_harris
