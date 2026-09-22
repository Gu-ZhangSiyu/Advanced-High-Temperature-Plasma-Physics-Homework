!***************************************************************
!   场数据工具子程序
!
!   调用示例
!   call Apply_Binomial_Filter(Jx, Jy, Jz, domain)
!   call Reduce_And_Normalize_Current(Jx, Jy, Jz, dx, dy, dz, domain)
!
!   子程序说明:
!   Apply_Binomial_Filter(Jx, Jy, Jz)
!       对电流场进行二项式平滑滤波，在YZ、XZ、XY三个平面分别做一维1-2-1平滑，
!       Jx, Jy, Jz				包含卫星层的电流数组 (0:nx+1, 0:ny+1, 0:nz+1)
!
!   Reduce_And_Normalize_Current(Jx, Jy, Jz, nx_in, ny_in, nz_in, dx_in, dy_in, dz_in)
!       累加相邻 rank 接口 ghost 的电流贡献，并除以网格体积归一化。
!       Jx, Jy, Jz				包含卫星层的电流数组 (0:nx+1, 0:ny+1, 0:nz+1)
!       nx_in, ny_in, nz_in		有效网格维度
!		dx_in, dy_in, dz_in		网格步长
!*********************************************************************
module Utility_Functions
    implicit none
contains
!*********************************************************************
subroutine Apply_Binomial_Filter(Jx, Jy, Jz, domain)
    use mesh_parameter, only: nx
    use Ghost_Cells_Module
    use Domain_Decomposition_Module, only: domain_decomposition
    implicit none
    real(8), intent(inout), contiguous :: Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    real(8), allocatable :: buf(:)
    integer :: i, k, max_dim

    max_dim = max(nx, domain%nz_local) + 2
    allocate(buf(0:max_dim))

    ! 2D XZ 中只沿 x、z 做平滑；y 单平面不再进行滤波。
    do k = 1, domain%nz_local; do i = 1, nx
        buf(0:domain%nz_local+1) = Jx(i,1,:)
        Jx(i,1,k) = 0.5d0*buf(k) + 0.25d0*(buf(k-1) + buf(k+1))
    end do; end do; call Ghost_Cells(Jx, 'Jx', domain)
    do k = 1, domain%nz_local; do i = 1, nx
        buf(0:nx+1) = Jx(:,1,k)
        Jx(i,1,k) = 0.5d0*buf(i) + 0.25d0*(buf(i-1) + buf(i+1))
    end do; end do; call Ghost_Cells(Jx, 'Jx', domain)

    do k = 1, domain%nz_local; do i = 1, nx
        buf(0:domain%nz_local+1) = Jy(i,1,:)
        Jy(i,1,k) = 0.5d0*buf(k) + 0.25d0*(buf(k-1) + buf(k+1))
    end do; end do; call Ghost_Cells(Jy, 'Jy', domain)
    do k = 1, domain%nz_local; do i = 1, nx
        buf(0:nx+1) = Jy(:,1,k)
        Jy(i,1,k) = 0.5d0*buf(i) + 0.25d0*(buf(i-1) + buf(i+1))
    end do; end do; call Ghost_Cells(Jy, 'Jy', domain)

    do k = 1, domain%nz_local; do i = 1, nx
        buf(0:domain%nz_local+1) = Jz(i,1,:)
        Jz(i,1,k) = 0.5d0*buf(k) + 0.25d0*(buf(k-1) + buf(k+1))
    end do; end do; call Ghost_Cells(Jz, 'Jz', domain)
    do k = 1, domain%nz_local; do i = 1, nx
        buf(0:nx+1) = Jz(:,1,k)
        Jz(i,1,k) = 0.5d0*buf(i) + 0.25d0*(buf(i-1) + buf(i+1))
    end do; end do; call Ghost_Cells(Jz, 'Jz', domain)
    deallocate(buf)
end subroutine Apply_Binomial_Filter
!*********************************************************************
subroutine Reduce_And_Normalize_Current(Jx, Jy, Jz, dx_in, dy_in, dz_in, domain)
	use mesh_parameter, only: nx, ny
    use Ghost_Cells_Module
    use Domain_Decomposition_Module, only: domain_decomposition, Accumulate_Current_Halos
    implicit none
    real(8), intent(in) :: dx_in, dy_in, dz_in
    real(8), intent(inout), contiguous :: Jx(0:,0:,0:), Jy(0:,0:,0:), Jz(0:,0:,0:)
    type(domain_decomposition), intent(in) :: domain
    real(8) :: node_vol_inv
    node_vol_inv = 1.0d0 / (dx_in * dy_in * dz_in)

    ! 接口ghost上可能含有邻居粒子的沉积贡献
    call Accumulate_Current_Halos(Jx, domain); call Accumulate_Current_Halos(Jy, domain); call Accumulate_Current_Halos(Jz, domain)
    Jx(1:nx, 1:ny, 1:domain%nz_local) = Jx(1:nx, 1:ny, 1:domain%nz_local) * node_vol_inv
    Jy(1:nx, 1:ny, 1:domain%nz_local) = Jy(1:nx, 1:ny, 1:domain%nz_local) * node_vol_inv
    Jz(1:nx, 1:ny, 1:domain%nz_local) = Jz(1:nx, 1:ny, 1:domain%nz_local) * node_vol_inv
    call Ghost_Cells(Jx, 'Jx', domain); call Ghost_Cells(Jy, 'Jy', domain); call Ghost_Cells(Jz, 'Jz', domain)
end subroutine Reduce_And_Normalize_Current
!*********************************************************************
end module Utility_Functions
