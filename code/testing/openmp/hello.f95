program parallel_hello_world
    use omp_lib

    !$omp parallel

    print *, 'hello from process: ', OMP_GET_THREAD_NUM()

    !$omp end parallel

end program
