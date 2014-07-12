#!/bin/bash

#cd $SAGE_ROOT

echo -e "\n$(tput setaf 4)Downloading the package SageManifolds v0.5...$(tput sgr 0)"
# downloading with either wget or curl 
wget -N http://users.camk.edu.pl/bejger/manifolds-0.5.tar.gz || curl -O http://users.camk.edu.pl/bejger/manifolds-0.5.tar.gz  

if [ ! -f "manifolds-0.5.tar.gz" ]; then
	echo  -e "\n$(tput setaf 1)Download the package manually and run again$(tput sgr 0)"
	exit 0 
fi 

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

# altering src/doc/en/reference/index.rst for documentation
if [[ -z $(grep "tensor_free_modules/index" src/doc/en/reference/index.rst) ]]; then
    
    echo "$(tput setaf 4)Adding tensor free modules documentation to index.rst...$(tput sgr 0)"
    sed -i '/<modules\/index>/a* :doc:`Tensors on free modules of finite rank <tensor_free_modules\/index>`' src/doc/en/reference/index.rst

fi

if [[ -z $(grep "manifolds/index" src/doc/en/reference/index.rst) ]]; then

    echo "$(tput setaf 4)Adding manifolds documentation to index.rst...$(tput sgr 0)"
    sed -i '/<tensor\/index>/a* :doc:`Differential Geometry <manifolds\/index>`' src/doc/en/reference/index.rst

fi

echo -e "\n$(tput bold)$(tput setaf 4)Running ./sage -b to re-build sage.$(tput sgr 0)"
./sage -b 

echo -e "\n$(tput bold)$(tput setaf 4)Building the html documentation (may take some time).$(tput sgr 0)"
./sage -docbuild reference inventory
./sage -docbuild reference html
./sage -docbuild tensors_free_module html
./sage -docbuild sagemanifolds html

echo echo -e "\n$(tput bold)$(tput setaf 4)Done!$(tput sgr 0)"

exit 0

