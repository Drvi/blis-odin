cd blis
make clean
make -j
cd ..
odin run example.odin -file -vet -strict-style -sanitize:address -debug