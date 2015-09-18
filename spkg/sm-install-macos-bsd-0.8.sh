#!/bin/bash

#cd $SAGE_ROOT

echo -e "\n$(tput setaf 4)Downloading the package SageManifolds 0.8...$(tput sgr 0)"
# Downloading with either wget or curl
wget -N http://sagemanifolds.obspm.fr/spkg/manifolds-0.8.tar.gz || curl -O http://sagemanifolds.obspm.fr/spkg/manifolds-0.8.tar.gz

if [ ! -f "manifolds-0.8.tar.gz" ]; then
	echo  -e "\n$(tput setaf 1)Download the package manually and run again$(tput sgr 0)"
	exit 0
fi

# Cleaning previous SM version if any
if [ -d  src/doc/en/reference/manifolds ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/doc/en/reference/manifolds$(tput sgr 0)"
    rm -fr src/doc/en/reference/manifolds/*
fi
if [ -d  src/doc/en/reference/tensor_free_modules ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/doc/en/reference/tensor_free_modules$(tput sgr 0)"
    rm -fr src/doc/en/reference/tensor_free_modules/*
fi
if [ -d src/sage/geometry/manifolds ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/sage/geometry/manifolds$(tput sgr 0)"
    rm -fr src/sage/geometry/manifolds/*
fi
if [ -d src/sage/tensor/modules ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/sage/tensor/modules$(tput sgr 0)"
    rm -fr src/sage/tensor/modules/*
fi

# Untaring the SM tree (src/sage/tensor/modules, src/sage/geometry/manifolds
# and the documentation files)
echo "$(tput setaf 4)Unpacking...$(tput sgr 0)"
tar -zxvf manifolds-0.8.tar.gz

# Altering src/sage/all.py for import
if [[ -z $(grep "tensor.modules.all" src/sage/all.py) ]]; then

   echo "$(tput setaf 4)Adding tensor.modules import call to all.py...$(tput sgr 0)"
   echo "from sage.tensor.modules.all import *" >> src/sage/all.py
fi

if [[ -z $(grep "geometry.manifolds.all" src/sage/all.py) ]]; then

   echo "$(tput setaf 4)Adding geometry.manifolds import call to all.py...$(tput sgr 0)"
   echo "from sage.geometry.manifolds.all import *" >> src/sage/all.py
fi

# Altering src/doc/en/reference/index.rst for documentation
if [[ -z $(grep "tensor_free_modules/index" src/doc/en/reference/index.rst) ]]; then

    echo "$(tput setaf 4)Adding tensor free modules documentation to index.rst...$(tput sgr 0)"
    gsed -i '/<modules\/index>/a* :doc:`Tensors on free modules of finite rank <tensor_free_modules\/index>`' src/doc/en/reference/index.rst

fi

if [[ -z $(grep "manifolds/index" src/doc/en/reference/index.rst) ]]; then

    echo "$(tput setaf 4)Adding manifolds documentation to index.rst...$(tput sgr 0)"
    gsed -i '/<tensor\/index>/a* :doc:`Differential Geometry <manifolds\/index>`' src/doc/en/reference/index.rst

fi

#
# Re-building sage
#
touch src/sage/tensor/modules/*.py
echo -e "\n$(tput bold)$(tput setaf 4)Running ./sage -b to re-build Sage.$(tput sgr 0)"
./sage -b

echo -e "\n$(tput bold)$(tput setaf 4)Building the html documentation (may take some time).$(tput sgr 0)"
./sage -docbuild reference inventory
./sage -docbuild reference html
echo -e "\n$(tput setaf 4)NB: if some error in building the documention occured,$(tput sgr 0)"
echo -e   "$(tput setaf 4)    (this may be caused by a previous installation of SageManifolds)$(tput sgr 0)"
echo -e   "$(tput setaf 4)    run make doc-clean && make doc$(tput sgr 0)"

echo -e "\n$(tput bold)$(tput setaf 4)Installation of SageManifolds 0.8 completed!\n$(tput sgr 0)"

exit 0
