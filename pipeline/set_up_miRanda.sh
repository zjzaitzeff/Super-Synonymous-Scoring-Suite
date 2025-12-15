cd miRanda/miRanda-3.3a

./configure
make CFLAGS="-fcommon"
cp src/miranda ../
cd ../..