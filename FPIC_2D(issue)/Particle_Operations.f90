!*********************************************************************
! 2D XZ 粒子推进与电流沉积。
!
! 空间采用 XZ 平面，电磁场仍为 2D3V：Ex/Ey/Ez、Bx/By/Bz 和
! Jx/Jy/Jz 均保留，粒子的 y 坐标固定为 y_plane，但 vy 参与 Boris
! 旋转并作为出平面电流 Jy 沉积。
!*********************************************************************
module Particle_Operations
    use mesh_parameter
    use pic_types
    use control_parameter, only: inv_dt
    use Domain_Decomposition_Module, only: domain_decomposition
    implicit none
contains
!*********************************************************************
subroutine Run_Particle_Pusher(part_array, nP, Ex, Ey, Ez, Bx, By, Bz, dt, rel_switch, qm_ratio_species, domain)
    use particle_solver_yee, only: Push_Particle
    implicit none
    type(particle), intent(inout) :: part_array(:)
    real(8), intent(in) :: Ex(0:, 0:, 0:), Ey(0:, 0:, 0:), Ez(0:, 0:, 0:)
    real(8), intent(in) :: Bx(0:, 0:, 0:), By(0:, 0:, 0:), Bz(0:, 0:, 0:)
    real(8), intent(in) :: dt, qm_ratio_species
    real(8) :: p_in(3), v_in(3), p_out(3), v_out(3), origin_vec(3), step_vec(3)
    integer, intent(in) :: nP, rel_switch
    type(domain_decomposition), intent(in) :: domain
    integer :: i, status_code
    integer :: dim_vec(3)

    origin_vec = [x_min, y_plane, domain%z_origin]
    step_vec = [dx, dy, dz]
    dim_vec = [nx, ny, domain%nz_local]
    do i = 1, nP
        if (part_array(i)%weight <= 0.0d0) cycle
        p_in = [part_array(i)%x, y_plane, part_array(i)%z]
        v_in = [part_array(i)%vx, part_array(i)%vy, part_array(i)%vz]

        ! 2D 物理边界只检查 x、z；y 不再是独立的粒子坐标。
        if (p_in(1) < x_min .or. p_in(1) >= x_max .or. &
            p_in(3) < z_min .or. p_in(3) >= z_max) then
            part_array(i)%weight = 0.0d0
            cycle
        end if

        call Push_Particle(pos_in=p_in, vel_in=v_in, qm=qm_ratio_species, dt=dt, &
                           Ex=Ex, Ey=Ey, Ez=Ez, Bx=Bx, By=By, Bz=Bz, &
                           grid_origin=origin_vec, grid_steps=step_vec, grid_dims=dim_vec, &
                           pos_out=p_out, vel_out=v_out, status=status_code, &
                           boundary_control=0, relativistic_control=rel_switch)
        if (status_code /= 0) then
            part_array(i)%weight = 0.0d0
        else
            part_array(i)%x = p_out(1); part_array(i)%y = y_plane; part_array(i)%z = p_out(3)
            part_array(i)%vx = v_out(1); part_array(i)%vy = v_out(2); part_array(i)%vz = v_out(3)
        end if
    end do
end subroutine Run_Particle_Pusher
!*********************************************************************
subroutine Deposit_Species_Yee(part, nP, Jx_out, Jy_out, Jz_out, species_charge, dt, domain, owned_cells_only)
    implicit none
    type(particle), intent(in) :: part(:)
    integer, intent(in) :: nP
    real(8), intent(inout) :: Jx_out(0:, 0:, 0:), Jy_out(0:, 0:, 0:), Jz_out(0:, 0:, 0:)
    real(8), intent(in) :: species_charge, dt
    type(domain_decomposition), intent(in) :: domain
    logical, optional, intent(in) :: owned_cells_only
    logical :: routed, intersects
    integer :: ip, iter, iu, iw, max_segments
    real(8) :: tu, tw, t_min, u1, w1, u2, w2, q_w_invdt, vy_norm
    real(8) :: u_curr, w_curr, u_next, w_next, u_min, u_max, w_min, w_max
    real(8), parameter :: eps = 1.0d-12

    routed = .false.
    if (present(owned_cells_only)) routed = owned_cells_only
    do ip = 1, nP
        if (part(ip)%weight <= 0.0d0) cycle
        q_w_invdt = species_charge * part(ip)%weight * inv_dt
        if (routed) q_w_invdt = species_charge * part(ip)%weight / dt

        ! 转换到 XZ 网格归一化坐标；y 位置固定，但 vy 仍产生 Jy。
        u2 = (part(ip)%x - x_min) * inv_dx
        w2 = (part(ip)%z - z_min) * inv_dz
        u1 = u2 - part(ip)%vx * dt * inv_dx
        w1 = w2 - part(ip)%vz * dt * inv_dz
        vy_norm = part(ip)%vy * dt * inv_dy

        if (routed) then
            ! 将一条长轨道裁剪到 physical box，保持完整宏观 dt 的电流。
            call Clip_Orbit(u1, w1, u2, w2, intersects)
            if (.not. intersects) cycle
        else
            if (u2 < 0.0d0 .or. u2 >= dble(nx-1) .or. &
                w2 < 0.0d0 .or. w2 >= dble(nz-1)) cycle
        end if
        max_segments = 4 + ceiling(abs(u2-u1)) + ceiling(abs(w2-w1))
        u_curr = u1; w_curr = w1

        do iter = 1, max_segments
            iu = floor(u_curr + eps)
            iw = floor(w_curr + eps)
            if (abs(u_curr - dble(iu)) < eps .and. u2 < u_curr) iu = iu - 1
            if (abs(w_curr - dble(iw)) < eps .and. w2 < w_curr) iw = iw - 1
            u_min = dble(iu); u_max = dble(iu + 1)
            w_min = dble(iw); w_max = dble(iw + 1)

            tu = 1.0d0; tw = 1.0d0
            if (u2 > u_max) tu = (u_max - u_curr) / (u2 - u_curr)
            if (u2 < u_min) tu = (u_min - u_curr) / (u2 - u_curr)
            if (w2 > w_max) tw = (w_max - w_curr) / (w2 - w_curr)
            if (w2 < w_min) tw = (w_min - w_curr) / (w2 - w_curr)
            t_min = min(tu, tw, 1.0d0)

            u_next = u_curr + t_min * (u2 - u_curr)
            w_next = w_curr + t_min * (w2 - w_curr)
            if (.not. routed .or. (iw+1 >= domain%k_start .and. iw+1 <= domain%k_end)) then
                call Deposit_Segment_VB(u_curr, w_curr, u_next, w_next, vy_norm, iu, iw, q_w_invdt, &
                                        Jx_out, Jy_out, Jz_out, domain)
            end if
            if (t_min >= 1.0d0) exit
            u_curr = u_next; w_curr = w_next
        end do
    end do
end subroutine Deposit_Species_Yee
!*********************************************************************
! 2D XZ Villasenor-Buneman 段沉积。
subroutine Deposit_Segment_VB(uA, wA, uB, wB, vy_norm, iu, iw, q_w_invdt, Jx_out, Jy_out, Jz_out, domain)
    use mesh_parameter, only: dx, dy, dz, nx, nz
    implicit none
    real(8), intent(in) :: uA, wA, uB, wB, vy_norm, q_w_invdt
    integer, intent(in) :: iu, iw
    integer :: i0, k0
    real(8), intent(inout) :: Jx_out(0:, 0:, 0:), Jy_out(0:, 0:, 0:), Jz_out(0:, 0:, 0:)
    type(domain_decomposition), intent(in) :: domain
    real(8) :: du, dw, u_mid, w_mid, Wx0, Wx1, Wz0, Wz1
    real(8) :: coef_x, coef_y, coef_z, Fx0, Fx1, Fy00, Fy01, Fy10, Fy11, Fz0, Fz1, du_dw_12

    du = uB - uA; dw = wB - wA
    du_dw_12 = (du * dw) / 12.0d0
    coef_x = q_w_invdt * dx; coef_y = q_w_invdt * dy; coef_z = q_w_invdt * dz
    u_mid = 0.5d0 * (uA + uB) - dble(iu)
    w_mid = 0.5d0 * (wA + wB) - dble(iw)
    i0 = iu + 1
    ! 全局 z 节点转为本地索引；接口贡献落到 ghost 平面并由邻居累加。
    k0 = (iw + 1) - domain%k_start + 1
    if (i0 < 1 .or. i0 >= nx) return
    if (iw + 1 < 1 .or. iw + 1 >= nz) return
    if (k0 < 0 .or. k0 > domain%nz_local) return

    Wx1 = u_mid; Wx0 = 1.0d0 - Wx1
    Wz1 = w_mid; Wz0 = 1.0d0 - Wz1
    ! x/z 面内电流只依赖二维轨迹；vy 的出平面电流用四节点权重。
    Fx0 = du * Wz0; Fx1 = du * Wz1
    Fy00 = vy_norm * (Wx0 * Wz0 + du_dw_12)
    Fy01 = vy_norm * (Wx0 * Wz1 - du_dw_12)
    Fy10 = vy_norm * (Wx1 * Wz0 - du_dw_12)
    Fy11 = vy_norm * (Wx1 * Wz1 + du_dw_12)
    Fz0 = dw * Wx0; Fz1 = dw * Wx1

    Jx_out(i0,   1, k0)   = Jx_out(i0,   1, k0)   + coef_x * Fx0
    Jx_out(i0,   1, k0+1) = Jx_out(i0,   1, k0+1) + coef_x * Fx1
    Jy_out(i0,   1, k0)   = Jy_out(i0,   1, k0)   + coef_y * Fy00
    Jy_out(i0,   1, k0+1) = Jy_out(i0,   1, k0+1) + coef_y * Fy01
    Jy_out(i0+1, 1, k0)   = Jy_out(i0+1, 1, k0)   + coef_y * Fy10
    Jy_out(i0+1, 1, k0+1) = Jy_out(i0+1, 1, k0+1) + coef_y * Fy11
    Jz_out(i0,   1, k0)   = Jz_out(i0,   1, k0)   + coef_z * Fz0
    Jz_out(i0+1, 1, k0)   = Jz_out(i0+1, 1, k0)   + coef_z * Fz1
end subroutine Deposit_Segment_VB
!*********************************************************************
! 按参数对 XZ 轨道进行裁剪。
subroutine Clip_Orbit(u1, w1, u2, w2, intersects)
    real(8), intent(inout) :: u1, w1, u2, w2
    logical, intent(out) :: intersects
    real(8) :: a(2), b(2), d(2), upper(2), t0, t1, ta, tb
    integer :: j

    a = [u1, w1]; b = [u2, w2]; d = b - a
    upper = [dble(nx-1), dble(nz-1)]
    t0 = 0.0d0; t1 = 1.0d0; intersects = .false.
    do j = 1, 2
        if (d(j) == 0.0d0) then
            if (a(j) < 0.0d0 .or. a(j) >= upper(j)) return
        else
            ta = -a(j) / d(j); tb = (upper(j) - a(j)) / d(j)
            t0 = max(t0, min(ta, tb)); t1 = min(t1, max(ta, tb))
        end if
    end do
    if (t1 <= t0) return
    b = a + t1*d; a = a + t0*d
    u1 = a(1); w1 = a(2); u2 = b(1); w2 = b(2)
    intersects = .true.
end subroutine Clip_Orbit
!*********************************************************************
! 向每个相交的 z-slab 发送临时轨迹记录。
subroutine Deposit_Long_Orbits(part, nP, Jx_out, Jy_out, Jz_out, species_charge, dt_macro, domain)
    use mpi_f08
    use Particle_Migration_Module, only: Pack_Particle, Unpack_Particle, packed_particle_size
    implicit none
    type(particle), intent(in) :: part(:)
    integer, intent(in) :: nP
    real(8), intent(inout) :: Jx_out(0:,0:,0:), Jy_out(0:,0:,0:), Jz_out(0:,0:,0:)
    real(8), intent(in) :: species_charge, dt_macro
    type(domain_decomposition), intent(in) :: domain
    integer :: sc(0:domain%nprocs-1), rc(0:domain%nprocs-1)
    integer :: sd(0:domain%nprocs-1), rd(0:domain%nprocs-1), cursor(0:domain%nprocs-1)
    real(8), allocatable :: sendbuf(:), recvbuf(:)
    type(particle) :: one(1)
    integer :: i, r, pass, lo, hi, slot, ierr
    real(8) :: za, zb

    sc = 0; sd = 0; rd = 0
    do pass = 1, 2
        do i = 1, nP
            if (part(i)%weight <= 0.0d0) cycle
            zb = (part(i)%z-z_min)*inv_dz
            za = zb-part(i)%vz*dt_macro*inv_dz
            lo = max(1, min(nz-1, floor(max(0.0d0, min(dble(nz-1), min(za,zb))))+1))
            hi = max(1, min(nz-1, floor(max(0.0d0, min(dble(nz-1), max(za,zb))))+1))
            do r = 0, domain%nprocs-1
                if (hi < domain%cuts(r) .or. lo >= domain%cuts(r+1)) cycle
                if (pass == 1) then
                    sc(r) = sc(r) + packed_particle_size
                else
                    slot = cursor(r) + 1
                    call Pack_Particle(part(i), sendbuf(slot:slot+packed_particle_size-1))
                    cursor(r) = cursor(r) + packed_particle_size
                end if
            end do
        end do
        if (pass == 1) then
            call MPI_Alltoall(sc, 1, MPI_INTEGER, rc, 1, MPI_INTEGER, domain%comm, ierr)
            do r = 1, domain%nprocs-1
                sd(r) = sd(r-1) + sc(r-1); rd(r) = rd(r-1) + rc(r-1)
            end do
            allocate(sendbuf(max(1,sum(sc))), recvbuf(max(1,sum(rc))))
            cursor = sd
        end if
    end do
    call MPI_Alltoallv(sendbuf, sc, sd, MPI_DOUBLE_PRECISION, recvbuf, rc, rd, &
                       MPI_DOUBLE_PRECISION, domain%comm, ierr)
    do slot = 1, sum(rc), packed_particle_size
        call Unpack_Particle(recvbuf(slot:slot+packed_particle_size-1), one(1))
        call Deposit_Species_Yee(one, 1, Jx_out, Jy_out, Jz_out, species_charge, dt_macro, domain, .true.)
    end do
    deallocate(sendbuf, recvbuf)
end subroutine Deposit_Long_Orbits
!*********************************************************************
end module Particle_Operations
