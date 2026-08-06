##########################################################
## ccAFv2: ccAFv2_Convert_onnx.py                       ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor,         ##
##           Kara Ramirez, and Thurston Herricks        ##
##                                                      ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

'''
Converting pytroch network into c code for R packages

This file is a example of how to convert a pytorch network into an onnx file.
The general steps are as follows.  Define the pytorch model so that flow control
is compatible with onnx.  This generally means that there are no if statements
or iterated lists in the forward function of the module.

The model below contains both the perceptron model and the calibrator layers to
scale the perceptron outputs to probabilites.


'''

import torch
import torch.nn as nn
from torch.fx.experimental.optimization import remove_dropout
import pathlib
import numpy as np

## Define ccAFv2 model

class ccAFv2_Static(nn.Module):
    def __init__(self):
        super().__init__()

        self.eps = torch.tensor(1e-10)

        # create network with activations.
        # Replicate ccAFv2 with nn.Sequential
        self.classifier = nn.Sequential(
            nn.Linear(in_features=871, out_features=871),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(in_features=871, out_features=871),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(in_features=871, out_features=7)
        )

        # create the calibrator
        self.calibrator = nn.Sequential(
            nn.Linear(in_features = 7, out_features = 7)
        )

        self._initialize_weights()


    # Add random weights to linear layers.
    def _initialize_weights(self):

        parent_path = pathlib.Path(__file__).parent

        # Load the weight dicts here.
        class_weight_path = parent_path / 'ccAFv2_RTorch_19.pth'
        #class_weights = torch.load(class_weight_path, weights_only = True)
        self.classifier.load_state_dict(dict(torch.load(class_weight_path, weights_only = True)))
        self.classifier = remove_dropout(self.classifier)
        self.classifier.eval()


        calib_weight_path = parent_path / 'ccAFv2_Rcorr_19.pth'
        #calib_weights = torch.load(calib_weight_path, weights_only = True)
        tmp_dict = dict(torch.load(calib_weight_path, weights_only = True))
        tmp_dict = {k.replace('network','0'):v for k,v in tmp_dict.items()}
        self.calibrator.load_state_dict(tmp_dict)
        self.calibrator.eval()

        '''
        #
        # add the weights to the classifier layers one at a time.  Since there
        # are dropouts layers during training, the linear layer is every third
        # layer in the class_weights dictionary.  The dropout layers are not
        # needed in this model.
        #
        self.classifier[0].weight = nn.Parameter(class_weights['0.weight'])
        self.classifier[0].bias   = nn.Parameter(class_weights['0.bias'])

        self.classifier[2].weight = nn.Parameter(class_weights['3.weight'])
        self.classifier[2].bias   = nn.Parameter(class_weights['3.bias'])

        self.classifier[4].weight = nn.Parameter(class_weights['6.weight'])
        self.classifier[4].bias   = nn.Parameter(class_weights['6.bias'])

        self.calibrator[0].weight  = nn.Parameter(calib_weights['network.weight'])
        self.calibrator[0].bias    = nn.Parameter(calib_weights['network.bias'])
        '''

    def forward(self, x):
        # calculate the logits from the input gene expression data
        x = self.classifier(x)
        # take the softmax of the logits
        x = torch.nn.functional.softmax(x, dim = 1)
        # ensure we don't have zeros when taking the log of the logits
        x = torch.log(x + self.eps)
        # apply the calibrator layer
        x = self.calibrator(x)
        # take the softmax of the calibrator layer and return the probabilities
        x = torch.nn.functional.softmax(x, dim = -1)

        return x


if __name__ == "__main__":

    # Get the path of the current folder
    parent_path = pathlib.Path(__file__).parent

    # assign the onnx output file name.
    onnx_path = parent_path / 'ccAFv2_r2.onnx'

    # create an input tensor to initialize the model.  We will process each s
    # sample individually without any batching in the R package.  This is slow
    # but saves on memory.
    input_tensor = torch.rand((1, 871), dtype=torch.float32)

    # Create the model, set it to eval, and run it once.
    ccAFv2_classifier = ccAFv2_Static()
    ccAFv2_classifier.eval()
    test = ccAFv2_classifier(input_tensor)

    # Export the onnx file
    torch.onnx.export(
        ccAFv2_classifier,      # model to export
        (input_tensor,),        # inputs of the model,
        onnx_path,              # filename of the ONNX model
        external_data = False,  # False onnx2c requires all data be in the onnx
        input_names=["input_tensor"],  # Rename inputs for the ONNX model
        dynamo=True             # True or False to select the exporter to use
    )



'''
Next steps to generate the c file.


docker run -it -v '/media/old_home/home/cplaisier/ccAFv2_test:/files' cplaisier/pytorch

pip install --upgrade torch

apt update
apt install git libprotobuf-dev protobuf-compiler

download and build the following repo:
https://github.com/kraiskil/onnx2c

git clone https://github.com/kraiskil/onnx2c.git
cd onnx2c
git submodule update --init

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make onnx2c

Then run the command:
(path to onnx2c)/onnx2c (path to dir)/ccAFv2_r2.onnx > (path to dir)/ccAFv2_r2.c

This will generate the c-code file of the network.  Please note there are probably
practical limits to the size of the network we can create.

The next steps are to add the file to the R package, compile the package and
test.

examples for adding c-code to a R package are here:
https://jbhender.github.io/Stats506/F17/Using_C_Cpp.html
https://jcarroll.com.au/2023/08/11/wrapping-c-code-in-an-r-package/


'''
