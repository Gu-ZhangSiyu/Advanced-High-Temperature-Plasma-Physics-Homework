!*********************************************************************
! 加权、有界、连续的z-slab重划分
!*********************************************************************
module Load_Balance_Module
    use mpi_f08
    use mesh_parameter, only: nx, ny, nz, z_min, inv_dz
    use pic_types, only: particle
    use control_parameter, only: sub_ratio, load_balance_cell_cost, load_balance_electron_cost, &
        load_balance_ion_cost, load_balance_threshold, load_balance_min_gain, load_balance_max_slab_factor
    use Domain_Decomposition_Module, only: domain_decomposition, Set_Domain_Cuts
    implicit none
    private
    public :: Propose_Balanced_Domain, Optimal_Slab_Cuts
contains
!*********************************************************************
! 最小-最大划分, 包含至少两个平面, 每个秩最多包含max_width个平面
subroutine Optimal_Slab_Cuts(cost, nprocs, max_width, cuts, peak)
    real(8), intent(in) :: cost(:)
    integer, intent(in) :: nprocs, max_width
    integer, intent(out) :: cuts(0:nprocs)
    real(8), intent(out) :: peak
    real(8) :: prefix(0:size(cost)), dp(0:size(cost)), next_dp(0:size(cost)), trial
    integer :: parent(nprocs,0:size(cost)), r, k, j, n, low, high
    n=size(cost)
    prefix(0)=0.0d0
    do k=1,n; prefix(k)=prefix(k-1)+cost(k); end do
    dp=huge(1.0d0); dp(0)=0.0d0; parent=-1
    do r=1,nprocs
        next_dp=huge(1.0d0)
        do k=2*r,min(r*max_width,n-2*(nprocs-r))
            low=max(2*(r-1),k-max_width); high=min((r-1)*max_width,k-2)
            do j=low,high
                trial=max(dp(j),prefix(k)-prefix(j))
                if (trial < next_dp(k)) then
                    next_dp(k)=trial; parent(r,k)=j
                end if
            end do
        end do; dp=next_dp
    end do
    peak=dp(n); cuts(nprocs)=n+1; k=n
    do r=nprocs,1,-1
        j=parent(r,k)
        if (j < 0) then
            cuts=-1
			return
        end if
        cuts(r-1)=j+1; k=j
    end do
end subroutine Optimal_Slab_Cuts
!*********************************************************************
subroutine Propose_Balanced_Domain(ions, ni, electrons, ne, old_domain, new_domain, changed, step)
    type(particle), intent(in) :: ions(:), electrons(:)
    integer, intent(in) :: ni, ne, step
    type(domain_decomposition), intent(in) :: old_domain
    type(domain_decomposition), intent(out) :: new_domain
    logical, intent(out) :: changed
    real(8) :: local_cost(nz), cost(nz), old_peak, new_peak, mean_cost
    integer :: i,k,r,ierr,max_width,cuts(0:old_domain%nprocs)
    new_domain=old_domain
    changed=.false.
    if (old_domain%nprocs == 1) return
    local_cost=0.0d0
    do i=1,ni
        if (ions(i)%weight <= 0.0d0) cycle
        k=max(1,min(nz-1,floor((ions(i)%z-z_min)*inv_dz)+1))
        local_cost(k)=local_cost(k)+load_balance_ion_cost/dble(sub_ratio)
    end do
    do i=1,ne
        if (electrons(i)%weight <= 0.0d0) cycle
        k=max(1,min(nz-1,floor((electrons(i)%z-z_min)*inv_dz)+1))
        local_cost(k)=local_cost(k)+load_balance_electron_cost
    end do
    call MPI_Allreduce(local_cost,cost,nz,MPI_DOUBLE_PRECISION,MPI_SUM,old_domain%comm,ierr)
    cost=cost+load_balance_cell_cost*dble(nx)*dble(ny)
    mean_cost=sum(cost)/dble(old_domain%nprocs)
    old_peak=0.0d0
    do r=0,old_domain%nprocs-1
        old_peak=max(old_peak,sum(cost(old_domain%cuts(r):old_domain%cuts(r+1)-1)))
    end do
    if (old_peak/mean_cost < load_balance_threshold) return
    max_width=min(nz,ceiling(load_balance_max_slab_factor* &
                  dble((nz+old_domain%nprocs-1)/old_domain%nprocs)))
    ! Rank 0 chooses cuts, 广播操作避免分区不一致
    if (old_domain%rank == 0) call Optimal_Slab_Cuts(cost,old_domain%nprocs,max_width,cuts,new_peak)
    call MPI_Bcast(cuts,size(cuts),MPI_INTEGER,0,old_domain%comm,ierr)
    call MPI_Bcast(new_peak,1,MPI_DOUBLE_PRECISION,0,old_domain%comm,ierr)
    if (any(cuts < 1)) call MPI_Abort(old_domain%comm,61,ierr)
    changed=any(cuts /= old_domain%cuts) .and. new_peak <= old_peak*(1.0d0-load_balance_min_gain)
    if (old_domain%rank == 0) then
        print '(A,I10,A,F9.4,A,F9.4,A,L1)', 'balance step=',step,' weighted max/avg=', &
            old_peak/mean_cost,' proposed=',new_peak/mean_cost,' accepted=',changed
    end if
    if (changed) call Set_Domain_Cuts(new_domain,cuts)
end subroutine Propose_Balanced_Domain
!*********************************************************************
end module Load_Balance_Module
