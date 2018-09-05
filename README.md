# SnapTree
fast tree data structure construction tool for big image data

# example

parallel convert data

```
time ./snaptree -i ./slice71 -o ./test71 -s 0 -e 32
time ./snaptree -i ./slice71 -o ./test71 -s 32 -e 64
time ./snaptree -i ./slice71 -o ./test71 -s 64 -e 71
```
generate mdata.bin only

```
time ./snaptree -i ./slice71 -o ./test71 -m 1 -x 42194 -y 54600 -z 71
```
