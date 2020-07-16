MPIWRAP_DEBUG=quiet \
#LD_PRELOAD=/opt/PARTIES/valgrind-3.12.0/lib/valgrind/libmpiwrap-amd64-linux.so \
mpirun -n 1 \
valgrind \
./parties > log.out 2> log.err &
