#!/bin/bash

#cd $SAGE_ROOT

echo -e "\n$(tput setaf 4)Downloading the package SageManifolds v0.5...$(tput sgr 0)"
wget -N http://sagemanifolds.obspm.fr/spkg/manifolds-0.5.tar.gz

# untaring the SM tree (tensor/modules, geometry/manifolds 
# and the documentation files)

echo "$(tput setaf 4)Unpacking...$(tput sgr 0)"
tar -zxvf manifolds-0.5.tar.gz

# altering the src/sage/all.py for import 

if [[ -z $(grep "tensor.modules.all" src/sage/all.py) ]]; then
   
   echo "$(tput setaf 4)Adding tensor.modules import call to all.py...$(tput sgr 0)"
   echo "from sage.tensor.modules.all import *" >> src/sage/all.py
fi

if [[ -z $(grep "geometry.manifolds.all" src/sage/all.py) ]]; then

   echo "$(tput setaf 4)Adding geometry.manifolds import call to all.py...$(tput sgr 0)"
   echo "from sage.geometry.manifolds.all import *" >> src/sage/all.py
fi 


echo -e "\n$(tput bold)$(tput setaf 4)Done! Run ./sage -b to re-build sage.$(tput sgr 0)"

exit 0

