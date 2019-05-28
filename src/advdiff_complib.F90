module advdiff_complib
  use advdiff_precision
  
  implicit none

  public :: sqrt_matrix, init_random_seed, randn
  public :: KCarte_to_KCanon, KCanon_to_KCarte
  public :: unfold_matrix
  
contains
  ! Random Sample from normal (Gaussian) distribution, Box-muller transform
  ! Source: https://sukhbinder.wordpress.com/fortran-random-number-generation/
   real(kind = dp) function randn()
    real(kind = dp) :: temp(2), r, theta
    real(kind = dp), parameter :: PI = 4.0_dp * atan (1.0_dp)
    
    call RANDOM_NUMBER(temp)
    r = (-2.0_dp*dlog(temp(1)))**0.5
    theta = 2.0_dp*PI*temp(2)
    randn = r*dsin(theta)
  end function randn
  
  subroutine init_random_seed()
    integer :: i, n
    integer, dimension(:), allocatable :: seed
    
    call RANDOM_SEED(size = n)
    allocate(seed(n))
          
    seed = 1989 + 611 + 37 * (/ (i - 1, i = 1, n) /)
    call RANDOM_SEED(PUT = seed)
          
    deallocate(seed)
  end subroutine
  
  ! Source: https://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
  subroutine sqrt_matrix(rootmat, mat)
    real(kind = dp), dimension(2,2), intent(out) :: rootmat
    real(kind = dp), dimension(2,2), intent(in) :: mat
    
    real(kind = dp) :: s, t, tau, delta
    
    tau = mat(1,1) + mat(2,2)
    delta = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
    
    s = dsqrt(delta)
    t = dsqrt(tau+2.0_dp*s)
    
    rootmat(1,1) = mat(1,1) + s
    rootmat(1,2) = mat(1,2)
    rootmat(2,1) = mat(2,1)
    rootmat(2,2) = mat(2,2) + s

    rootmat = rootmat/t
  end subroutine

  ! Directly translate from MATLAB code
  subroutine KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)
    real(kind = dp), intent(out) :: sigma1_sq, sigma2_sq, phi_K
    real(kind = dp), intent(in) :: Kxx, Kyy, Kxy
    
    real(kind = dp) :: sigma_sq_sum, sigma_sq_diff
    
    sigma_sq_sum = Kxx + Kyy
    sigma_sq_diff = dsqrt((2.0_dp*Kxy)**2 + (Kxx-Kyy)**2)
    
    phi_K = 0.5_dp*datan2(2.0_dp*Kxy, Kxx-Kyy)
    
    sigma1_sq = 0.5_dp*(sigma_sq_sum + sigma_sq_diff)
    sigma2_sq = 0.5_dp*(sigma_sq_sum - sigma_sq_diff)
  end subroutine

    ! Directly translate from MATLAB code
  subroutine KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)
    real(kind = dp), intent(in) :: sigma1_sq, sigma2_sq, phi_K
    real(kind = dp), intent(out) :: Kxx, Kyy, Kxy
    
    real(kind = dp) :: c, s, sigma_sq_sumh, sigma_sq_diffh
    
    c = dcos(2.0_dp*phi_K)
    s = dsin(2.0_dp*phi_K)

    sigma_sq_sumh = 0.5_dp*(sigma1_sq + sigma2_sq)
    sigma_sq_diffh = 0.5_dp*(sigma1_sq - sigma2_sq)

    Kxx = sigma_sq_sumh + sigma_sq_diffh*c
    Kyy = sigma_sq_sumh - sigma_sq_diffh*c
    Kxy = sigma_sq_diffh*s
  end subroutine
  
  subroutine unfold_matrix(vec, arr2d)
    real(kind=dp), dimension(:), intent(inout) :: vec
    real(kind=dp), dimension(:,:), intent(in) :: arr2d
    
    integer :: i, j, counter
    
    if (size(vec, 1) .ne. (size(arr2d,1)*size(arr2d, 2))) then
      write(6, *) "Inconsistent dimension of arr2d, vec"
      stop 1
    end if
    
    counter = 0
    do j = 1, size(arr2d, 2)
      do i = 1, size(arr2d, 1)
        counter = counter + 1
        vec(counter) = arr2d(i, j)
      end do
    end do

  end subroutine unfold_matrix
  
end module advdiff_complib
