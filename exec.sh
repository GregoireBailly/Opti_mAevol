#-DCMAKE_BUILD_TYPE=Release
cmake . -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=on
make
datestr=$(date +%Y-%m-%d_%H:%M:%S)
mkdir simulation_$datestr
cd simulation_$datestr
export OMP_DISPLAY_ENV=TRUE
../pdc_micro_aevol -n 500 -w 32 -h 32
