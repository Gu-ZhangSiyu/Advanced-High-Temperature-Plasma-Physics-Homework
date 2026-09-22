!*********************************************************************
!   填充三维数组的卫星层(ghost cells)
!
!   调用示例
!   call Ghost_Cells(F, 'Bx', domain)
!
!   子程序说明:
!   Ghost_Cells(F, var_name)
!       填充三维数组的卫星层(ghost cells)。
!       径向(X): 理想导电壁(PEC) 奇偶性边界
!       y方向  : 2D XZ 单平面复制，仅用于插值/存储保护，不代表物理边界
!       轴向(Z)  : rank 内部交换 halo，全局端点采用零梯度出流
!       输入/输出: F — 包含卫星层的三维数组
!       输入: var_name — 变量名字符串，决定边界类型
!*********************************************************************
module Ghost_Cells_Module
    use Domain_Decomposition_Module, only: domain_decomposition, Exchange_Field_Halos
    implicit none
contains
!*********************************************************************
subroutine Ghost_Cells(F, var_name, domain)
    use mesh_parameter, only: nx, ny
    implicit none
    real(8), intent(inout), contiguous :: F(0:, 0:, 0:)
    character(len=*), intent(in) :: var_name
    type(domain_decomposition), intent(in) :: domain
    select case (trim(var_name))
        case ('Ex', 'Jx')
            F(0, :, :) = F(1, :, :); F(nx+1, :, :) = F(nx, :, :)
        case ('Ey', 'Jy')
            F(0, :, :) = -F(1, :, :); F(nx+1, :, :) = -F(nx, :, :)
        case ('Ez', 'Jz')
            F(0, :, :) = -F(1, :, :); F(nx+1, :, :) = -F(nx, :, :)
        case ('Bx')
            F(0, :, :) = -F(1, :, :); F(nx+1, :, :) = -F(nx, :, :)
        case ('By')
            F(0, :, :) = F(1, :, :); F(nx+1, :, :) = F(nx, :, :)
        case ('Bz')
            F(0, :, :) = F(1, :, :); F(nx+1, :, :) = F(nx, :, :)
        case default
            print *, "ERROR: Invalid variable passed to Ghost_Cells"
            stop
    end select
    ! y 方向已经降为单平面；两侧只复制同一层，避免任何隐含的 y 导数。
    F(:, 0, :) = F(:, 1, :); F(:, ny+1, :) = F(:, ny, :)
    ! z向内部边界由邻居halo交换填充, 全局端点保持零梯度开放边界。
    call Exchange_Field_Halos(F, domain)
end subroutine Ghost_Cells
!*********************************************************************
end module Ghost_Cells_Module
