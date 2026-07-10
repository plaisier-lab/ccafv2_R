# Instructions to use onnx output from pytorch to compile C code for the model and calibrator

## Making the onnx file
Using a file based on the below convert the models into an onnx file for conversion:

### Compile the pytorch model in onnx format
In Bash:
```
docker run -it -v '/media/old_home/home/cplaisier/:/files' cplaisier/pytorch
pip install --upgrade torch
pip install onnxscript
cd <location of model to compile>
python ccAFv2_Convert_onnx.py
```

### Setup for making the C code
In Bash:
```
docker run -it -v '/media/old_home/home/cplaisier:/files' cplaisier/pytorch
pip install --upgrade torch
apt update
apt install git libprotobuf-dev protobuf-compiler
```

### Compile onnx2c
In Bash:
```
cd /opt
git clone https://github.com/kraiskil/onnx2c.git
cd onnx2c
git submodule update --init
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make onnx2c
```

### Making the C code
Change directory to the where the onnx file is output:
In Bash:
```
/opt/onnx2c/build/onnx2c ccAFv2_r2.onnx > C_ccAFv2_r2.c
```
