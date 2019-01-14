!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine unique_id_rng_seed(summ)
    use :: iso_c_binding
    integer :: seed_size
    integer(kind = 8), intent(out) :: summ
    integer, allocatable :: sseedd(:)

    call random_seed() ! initialize with system generated seed
    call random_seed(size = seed_size) ! find out size of seed
    allocate(sseedd(seed_size))
    call random_seed(get = sseedd) ! get system generated seeds
    !write(*,*) sseedd            ! writes system generated seeds

    summ = 0
    do n = 1, seed_size
        summ = summ + sseedd(n);
    enddo

end subroutine unique_id_rng_seed