gcc -O3 -Wall -lpthread -lm -g read_barcodelist.c fast_map.c unfold_map.c -o unfold_map
gcc -O3 -Wall shm_initialize.c -o shm_initialize
gcc -O3 -Wall shm_release.c -o shm_release
gcc -O3 -Wall -lm -g lowess.c stats.c -o unfold_stats
