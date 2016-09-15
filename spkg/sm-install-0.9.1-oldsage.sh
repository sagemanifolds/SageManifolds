#!/bin/bash

if [ "$1" = "doc" ]; then
    echo -e "\n$(tput bold)$(tput setaf 4)Building the html documentation (may take some time).$(tput sgr 0)"
    ./sage -docbuild reference/parallel html
    ./sage -docbuild reference/tensor_free_modules html
    ./sage -docbuild reference/manifolds html
    echo -e "\n$(tput setaf 4)NB: if some error in building the documention occured,$(tput sgr 0)"
    echo -e   "$(tput setaf 4)    (this may be caused by a previous installation of SageManifolds)$(tput sgr 0)"
    echo -e   "$(tput setaf 4)    run the following command:$(tput sgr 0)"
    echo -e   "$(tput setaf 4)        make doc-clean && make doc$(tput sgr 0)"
    exit 0
fi

if [ "$1" != "no-download" ]; then
    echo -e "\n$(tput setaf 4)Downloading the sources of SageManifolds 0.9.1...$(tput sgr 0)"
    # Downloading with either wget or curl
    wget -N http://sagemanifolds.obspm.fr/spkg/manifolds-0.9.1.tar.gz || curl -O http://sagemanifolds.obspm.fr/spkg/manifolds-0.9.1.tar.gz

    if [ ! -f "manifolds-0.9.1.tar.gz" ]; then
        echo  -e "\n$(tput setaf 1)Download the package manually and run again$(tput sgr 0)"
        exit 0
    fi
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
if [ -d src/sage/manifolds ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/sage/manifolds$(tput sgr 0)"
    rm -fr src/sage/manifolds/*
fi
if [ -d src/sage/geometry/manifolds ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/sage/geometry/manifolds$(tput sgr 0)"
    rm -fr src/sage/geometry/manifolds
fi
if [ -d src/sage/tensor/modules ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/sage/tensor/modules$(tput sgr 0)"
    rm -fr src/sage/tensor/modules/*
fi
if [[ -n $(grep "sage.geometry.manifolds.all" src/sage/all.py) ]]; then
   echo "$(tput setaf 4)Delete old manifolds import from all.py...$(tput sgr 0)"
   sed -i '/^from sage.geometry.manifolds.all/d' src/sage/all.py || \
   gsed -i '/^from sage.geometry.manifolds.all/d' src/sage/all.py
fi
if [[ -n $(grep "sage.tensor.modules.all" src/sage/all.py) ]]; then
   echo "$(tput setaf 4)Delete old tensor modules import from all.py...$(tput sgr 0)"
   sed -i '/^from sage.tensor.modules.all/d' src/sage/all.py || \
   gsed -i '/^from sage.tensor.modules.all/d' src/sage/all.py
fi

# Untaring the SM tree (src/sage/manifolds, src/sage/tensor/modules,
# src/sage/parallel and the documentation files)
echo "$(tput setaf 4)Unpacking...$(tput sgr 0)"
tar -zxvf manifolds-0.9.1.tar.gz

# Altering src/sage/all.py for import
if [[ -z $(grep "sage.manifolds.all" src/sage/all.py) ]]; then
   echo "$(tput setaf 4)Adding manifolds import to all.py...$(tput sgr 0)"
   echo "from sage.manifolds.all import *" >> src/sage/all.py
fi

# Altering src/doc/en/reference/index.rst for documentation
# using either sed (linux) or gsed (MacOS X)
if [[ -z $(grep "tensor_free_modules/index" src/doc/en/reference/index.rst) ]]; then
    echo "$(tput setaf 4)Adding tensor free modules documentation to index.rst...$(tput sgr 0)"
    sed -i '/<modules\/index>/a* :doc:`Tensors on free modules of finite rank <tensor_free_modules\/index>`' src/doc/en/reference/index.rst || \
    gsed -i '/<modules\/index>/a* :doc:`Tensors on free modules of finite rank <tensor_free_modules\/index>`' src/doc/en/reference/index.rst
fi

if [[ -z $(grep "manifolds/index" src/doc/en/reference/index.rst) ]]; then
    echo "$(tput setaf 4)Adding manifolds documentation to index.rst...$(tput sgr 0)"
    sed -i '/<tensor\/index>/a* :doc:`Manifolds <manifolds\/index>`' src/doc/en/reference/index.rst || \
    gsed -i '/<tensor\/index>/a* :doc:`Manifolds <manifolds\/index>`' src/doc/en/reference/index.rst
fi

#
# Re-building sage
#
touch src/sage/manifolds/*.py
touch src/sage/manifolds/differentiable/*.py
touch src/sage/tensor/modules/*.py
touch src/sage/parallel/*.py
echo -e "\n$(tput bold)$(tput setaf 4)Running ./sage -b to re-build Sage.$(tput sgr 0)"
./sage -b
echo -e "\n$(tput bold)$(tput setaf 4)Installation of SageManifolds 0.9.1 completed!\n$(tput sgr 0)"

exit 0
