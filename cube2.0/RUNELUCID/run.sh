source module_load_t7920.sh

cd ./utilities/
make clean
make 

cd ../
make clean
./utilities/ic.x

make && ./main.x && ./utilities/cicpower.x  
&& ./utilities/idsp.x
#&& ./utilities/fof.x
