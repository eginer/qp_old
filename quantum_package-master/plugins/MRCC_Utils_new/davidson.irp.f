subroutine davidson_diag_mrcc(dets_in,u_in,energies,dim_in,sze,N_st,Nint,iunit,istate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization.
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! iunit : Unit number for the I/O
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, Nint, iunit, istate
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision, intent(inout) :: u_in(dim_in,N_st)
  double precision, intent(out)  :: energies(N_st)
  double precision, allocatable  :: H_jj(:)
  
  double precision               :: diag_h_mat_elem
  integer                        :: i
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_bielec_integrals_in_map
  allocate(H_jj(sze))
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(sze,H_jj,N_det_ref,dets_in,Nint,istate,delta_ii,idx_ref)           &
      !$OMP  PRIVATE(i)
  !$OMP DO SCHEDULE(guided)
  do i=1,sze
    H_jj(i) = diag_h_mat_elem(dets_in(1,1,i),Nint) 
  enddo
  !$OMP END DO 
  !$OMP DO SCHEDULE(guided)
  do i=1,N_det_ref
    H_jj(idx_ref(i)) +=  delta_ii(i,istate)
  enddo
  !$OMP END DO 
  !$OMP END PARALLEL

  call davidson_diag_hjj_mrcc(dets_in,u_in,H_jj,energies,dim_in,sze,N_st,Nint,iunit,istate)
  deallocate (H_jj)
end

subroutine davidson_diag_hjj_mrcc(dets_in,u_in,H_jj,energies,dim_in,sze,N_st,Nint,iunit,istate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! iunit : Unit for the I/O
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, Nint, istate
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision,  intent(in)  :: H_jj(sze)
  integer,  intent(in)  :: iunit
  double precision, intent(inout) :: u_in(dim_in,N_st)
  double precision, intent(out)  :: energies(N_st)
  
  integer                        :: iter
  integer                        :: i,j,k,l,m
  logical                        :: converged
  
  double precision               :: overlap(N_st,N_st)
  double precision               :: u_dot_v, u_dot_u
  
  integer, allocatable           :: kl_pairs(:,:)
  integer                        :: k_pairs, kl
  
  integer                        :: iter2
  double precision, allocatable  :: W(:,:,:),  U(:,:,:), R(:,:)
  double precision, allocatable  :: y(:,:,:,:), h(:,:,:,:), lambda(:)
  double precision               :: diag_h_mat_elem
  double precision               :: residual_norm(N_st)
  character*(16384)              :: write_buffer
  double precision               :: to_print(2,N_st)
  double precision               :: cpu, wall
  

  call write_time(iunit)
  call wall_time(wall)
  call cpu_time(cpu)
  write(iunit,'(A)') ''
  write(iunit,'(A)') 'Davidson Diagonalization'
  write(iunit,'(A)') '------------------------'
  write(iunit,'(A)') ''
  call write_int(iunit,N_st,'Number of states')
  call write_int(iunit,sze,'Number of determinants')
  write(iunit,'(A)') ''
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================'
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write_buffer = ' Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'           Energy         Residual'
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================'
  enddo
  write(iunit,'(A)') trim(write_buffer)

  allocate(                                                          &
      kl_pairs(2,N_st*(N_st+1)/2),                                   &
      W(sze,N_st,davidson_sze_max),                                                   &
      U(sze,N_st,davidson_sze_max),                                  &
      R(sze,N_st),                                                   &
      h(N_st,davidson_sze_max,N_st,davidson_sze_max),                &
      y(N_st,davidson_sze_max,N_st,davidson_sze_max),                &
      lambda(N_st*davidson_sze_max))
  
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  
  ! Initialization
  ! ==============
  
  k_pairs=0
  do l=1,N_st
    do k=1,l
      k_pairs+=1
      kl_pairs(1,k_pairs) = k
      kl_pairs(2,k_pairs) = l
    enddo
  enddo
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(U,sze,N_st,overlap,kl_pairs,k_pairs,             &
      !$OMP  Nint,dets_in,u_in)                                 &
      !$OMP  PRIVATE(k,l,kl,i)
  
  
  ! Orthonormalize initial guess
  ! ============================
  
  !$OMP DO
  do kl=1,k_pairs
    k = kl_pairs(1,kl)
    l = kl_pairs(2,kl)
    if (k/=l) then
      overlap(k,l) = u_dot_v(U_in(1,k),U_in(1,l),sze)
      overlap(l,k) = overlap(k,l)
    else
      overlap(k,k) = u_dot_u(U_in(1,k),sze)
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call ortho_lowdin(overlap,size(overlap,1),N_st,U_in,size(U_in,1),sze)
  
  ! Davidson iterations
  ! ===================
  
  converged = .False.
  
  do while (.not.converged)
    
    !$OMP PARALLEL DEFAULT(NONE)                                     &
        !$OMP PRIVATE(k,i) SHARED(U,u_in,sze,N_st)
    do k=1,N_st
      !$OMP DO
      do i=1,sze
        U(i,k,1) = u_in(i,k)
      enddo
      !$OMP END DO 
    enddo
    !$OMP END PARALLEL
    
    do iter=1,davidson_sze_max-1
      
      ! Compute W_k = H |u_k>
      ! ----------------------
      
      do k=1,N_st
        call H_u_0_mrcc(W(1,k,iter),U(1,k,iter),H_jj,sze,dets_in,Nint,istate)
      enddo
      
      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      do l=1,N_st
        do k=1,N_st
          do iter2=1,iter-1
            h(k,iter2,l,iter) = u_dot_v(U(1,k,iter2),W(1,l,iter),sze)
            h(k,iter,l,iter2) = h(k,iter2,l,iter)
          enddo
        enddo
        do k=1,l
          h(k,iter,l,iter) = u_dot_v(U(1,k,iter),W(1,l,iter),sze)
          h(l,iter,k,iter) = h(k,iter,l,iter)
        enddo
      enddo

      !DEBUG H MATRIX
      !do i=1,iter
      !  print '(10(x,F16.10))',  h(1,i,1,1:i)
      !enddo
      !print *,  ''
      !END
      
      ! Diagonalize h
      ! -------------
      call lapack_diag(lambda,y,h,N_st*davidson_sze_max,N_st*iter)
      
      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------
      
      do k=1,N_st
        do i=1,sze
          U(i,k,iter+1) = 0.d0
          W(i,k,iter+1) = 0.d0
          do l=1,N_st
            do iter2=1,iter
              U(i,k,iter+1) = U(i,k,iter+1) + U(i,l,iter2)*y(l,iter2,k,1)
              W(i,k,iter+1) = W(i,k,iter+1) + W(i,l,iter2)*y(l,iter2,k,1)
            enddo
          enddo
        enddo
      enddo
      
      ! Compute residual vector
      ! -----------------------
      
      do k=1,N_st
        do i=1,sze
          R(i,k) = lambda(k) * U(i,k,iter+1) - W(i,k,iter+1)
        enddo
        residual_norm(k) = u_dot_u(R(1,k),sze)
        to_print(1,k) = lambda(k) + nuclear_repulsion
        to_print(2,k) = residual_norm(k)
      enddo

      write(iunit,'(X,I3,X,100(X,F16.10,X,E16.6))'), iter, to_print(:,1:N_st)
      call davidson_converged(lambda,residual_norm,wall,iter,cpu,N_st,converged)
      if (converged) then
        exit
      endif
      
      
      ! Davidson step
      ! -------------
      
      do k=1,N_st
        do i=1,sze
          U(i,k,iter+1) = -1.d0/max(H_jj(i) - lambda(k),1.d-2) * R(i,k)
        enddo
      enddo
      
      ! Gram-Schmidt
      ! ------------
      
      double precision               :: c
      do k=1,N_st
        do iter2=1,iter
          do l=1,N_st
            c = u_dot_v(U(1,k,iter+1),U(1,l,iter2),sze)
            do i=1,sze
              U(i,k,iter+1) -= c * U(i,l,iter2)
            enddo
          enddo
        enddo
        do l=1,k-1
          c = u_dot_v(U(1,k,iter+1),U(1,l,iter+1),sze)
          do i=1,sze
            U(i,k,iter+1) -= c * U(i,l,iter+1)
          enddo
        enddo
        call normalize( U(1,k,iter+1), sze )
      enddo
    
      !DEBUG : CHECK OVERLAP
      !print *,  '==='
      !do k=1,iter+1
      !  do l=1,k
      !  c = u_dot_v(U(1,1,k),U(1,1,l),sze)
      !  print *,  k,l, c
      !  enddo
      !enddo
      !print *,  '==='
      !pause
      !END DEBUG


    enddo

    if (.not.converged) then
      iter = davidson_sze_max-1
    endif
    
    ! Re-contract to u_in
    ! -----------
    
    do k=1,N_st
      energies(k) = lambda(k)
      do i=1,sze
        u_in(i,k) = 0.d0
        do iter2=1,iter
          do l=1,N_st
            u_in(i,k) += U(i,l,iter2)*y(l,iter2,k,1)
          enddo
        enddo
      enddo
    enddo

  enddo

  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================'
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write(iunit,'(A)') ''
  call write_time(iunit)

  deallocate (                                                       &
      kl_pairs,                                                      &
      W,                                                             &
      U,                                                             &
      R,                                                             &
      h,                                                             &
      y,                                                             &
      lambda                                                         &
      )
  abort_here = abort_all
end

subroutine H_u_0_mrcc(v_0,u_0,H_jj,n,keys_tmp,Nint,istate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0>
  !   
  ! n : number of determinants
  !     
  ! H_jj : array of <j|H|j>
  END_DOC
  integer, intent(in)            :: n,Nint,istate
  double precision, intent(out)  :: v_0(n)
  double precision, intent(in)   :: u_0(n)
  double precision, intent(in)   :: H_jj(n)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  integer, allocatable           :: idx(:)
  double precision               :: hij
  double precision, allocatable  :: vt(:)
  integer                        :: i,j,k,l, jj,ii
  integer                        :: i0, j0
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (n>0)
  PROVIDE ref_bitmask_energy delta_ij 
  integer, parameter             :: block_size = 157
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE(i,hij,j,k,idx,jj,ii,vt)                             &
      !$OMP SHARED(n_det_ref,n_det_non_ref,idx_ref,idx_non_ref,n,H_jj,u_0,keys_tmp,Nint,v_0,istate,delta_ij)
  !$OMP DO SCHEDULE(static)
  do i=1,n  
    v_0(i) = H_jj(i) * u_0(i)
  enddo 
  !$OMP END DO
  allocate(idx(0:n), vt(n))
  Vt = 0.d0
  !$OMP DO SCHEDULE(guided)
  do i=1,n
    idx(0) = i
    call filter_connected_davidson(keys_tmp,keys_tmp(1,1,i),Nint,i-1,idx)
    do jj=1,idx(0)
      j = idx(jj)
      if ( (dabs(u_0(j)) > 1.d-7).or.((dabs(u_0(i)) > 1.d-7)) ) then
        call i_H_j(keys_tmp(1,1,j),keys_tmp(1,1,i),Nint,hij)
        hij = hij 
        vt (i) = vt (i) + hij*u_0(j)
        vt (j) = vt (j) + hij*u_0(i)
      endif
    enddo
  enddo
  !$OMP END DO

  !$OMP DO SCHEDULE(guided)
  do ii=1,n_det_ref
    i = idx_ref(ii)
    do jj = 1, n_det_non_ref
        j = idx_non_ref(jj)
        vt (i) = vt (i) + delta_ij(ii,jj,istate)*u_0(j)
        vt (j) = vt (j) + delta_ij(ii,jj,istate)*u_0(i)
    enddo
  enddo
  !$OMP END DO
  !$OMP CRITICAL
  do i=1,n
    v_0(i) = v_0(i) + vt(i)
  enddo
  !$OMP END CRITICAL
  deallocate(idx,vt)
  !$OMP END PARALLEL
end


