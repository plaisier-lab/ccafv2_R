# Instructions to use onnx output from pytorch to compile C code for the model and calibrator

## Making the onnx file
Using a file based on the below convert the models into an onnx file for conversion:

### Step 1: Compile the pytorch model in onnx format
In Bash:
```bash
docker run -it -v '/media/old_home/home/cplaisier/:/files' cplaisier/pytorch
pip install --upgrade torch
pip install onnxscript
cd <location of model to compile>
python ccAFv2_Convert_onnx.py
```

### Step 2: Setup for making the C code
In Bash:
```bash
docker run -it -v '/media/old_home/home/cplaisier:/files' cplaisier/pytorch
pip install --upgrade torch
apt update
apt install git libprotobuf-dev protobuf-compiler
```

### Step 3: Compile onnx2c
In Bash:
```bash
cd /opt
git clone https://github.com/kraiskil/onnx2c.git
cd onnx2c
git submodule update --init
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make onnx2c
```

### Step 4: Making the C code
Change directory to the where the onnx file is output:
In Bash:
```bash
/opt/onnx2c/build/onnx2c ccAFv2_r2.onnx > C_ccAFv2_r2.c
```

### Step 5: Modify the C code to allow R embedding
Add this to the top of the C file with other include statements:
```c
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
```

Add this to the end of the C file for the R callout:
```c
// Wrapper for R in C function
SEXP C_ccAFv2(SEXP inp_vec) {
    if (!isReal(inp_vec) || LENGTH(inp_vec) != 871) {
        error("Input must be a numeric vector of length 871");
    }

    // Create pointer to the input vector
    const double* r_inp_ptr = REAL(inp_vec);
    // Create a float array for the input data and output data arrays
    float inp_c[1][871];
    float oup_c[1][7];

    // Convert input vector from double to float
    for (int i = 0; i < 871; i++) {
        inp_c[0][i] = (float) r_inp_ptr[i];
    }

    // Call classifier of vector
    entry(inp_c, oup_c);
    
    //Prepair the S-expression for R 
    SEXP result = PROTECT(Rf_allocVector(REALSXP, 7));
    double* result_ptr = REAL(result);
   
    // Copy the output data array to the S-expression
    for (int i = 0; i < 7; i++){
        result_ptr[i] = (double) oup_c[0][i];
    }

    UNPROTECT(1);
    return result;
}
```

### Step 6: Embed C code into R package

#### Part A: Modify R package NAMESPACE file
Modify NAMESPACE file by adding something like the following:
```R
export(ccAFv2_classifier)
useDynLib(ccAFv2, .registration = TRUE)
```
#### Part B: Make a src directory and add two c files
- File 1: 'C_ccAFv2.C' from above which is the model
- File 2: 'init.c' which should look like the below:
```c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP C_ccAFv2(SEXP inp_vec);

static const R_CallMethodDef CallEntries[] = {
    {"C_ccAFv2", (DL_FUNC) &C_ccAFv2, 1},
    {NULL, NULL, 0}
};

void R_init_ccAFv2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
```

#### Part C: Modify R package code to add in the function to call the classifier
Add this R code to call the compiled C code classifier:
```R
#' ccAFv2 classifier function
#' 
#' This function calls the C code C_ccAFv2.c to run the neural network classifier. 
#' ccAFv2_classifier returns a 1 x 7 double precision vector of predictions 
#' 
#' @param  norm_expVec a double precision 1 x 861 vector of normalized gene expression values
#' @return  ccAFv2_classifier returns a 1 x 7 double precision vector of class probabilites
#'
#' @export
ccAFv2_classifier = function(norm_expVec) {
    .Call("C_ccAFv2", norm_expVec)
}
```
Then the function ccAFv2_classifier will run the C code using the input and return an output result tensor/vector.

