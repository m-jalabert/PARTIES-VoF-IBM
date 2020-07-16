MPIWRAP_DEBUG=quiet \
LD_PRELOAD=/opt/valgrind-3.12.0/lib/valgrind/libmpiwrap-amd64-linux.so \
mpirun -n 8 \
valgrind --leak-check=full \
./parties > log.out 2> log.err &
