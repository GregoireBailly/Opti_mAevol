#-DCMAKE_BUILD_TYPE=Release
cmake . -DCMAKE_BUILD_TYPE=Release
make
datestr=$(date +%Y-%m-%d_%H:%M:%S)
mkdir simulation_$datestr
cd simulation_$datestr
../pdc_micro_aevol -n 50 -w 32 -h 32
