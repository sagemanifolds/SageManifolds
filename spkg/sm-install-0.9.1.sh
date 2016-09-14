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

# Cleaning previous SageManifolds version if any
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
if [ -d src/sage/tensor/modules ]; then
	echo  -e "$(tput setaf 4)Removing previous version of src/sage/tensor/modules$(tput sgr 0)"
    rm -fr src/sage/tensor/modules/*
fi

# Untaring the SageManifolds tree (src/sage/manifolds, src/sage/tensor/modules,
# src/sage/parallel and the documentation files)
echo "$(tput setaf 4)Unpacking...$(tput sgr 0)"
tar -zxvf manifolds-0.9.1.tar.gz

#
# Re-building Sage
#
touch src/sage/manifolds/*.py
touch src/sage/manifolds/differentiable/*.py
touch src/sage/tensor/modules/*.py
touch src/sage/parallel/*.py
echo -e "\n$(tput bold)$(tput setaf 4)Running ./sage -b to re-build Sage.$(tput sgr 0)"
./sage -b
echo -e "\n$(tput bold)$(tput setaf 4)Installation of SageManifolds 0.9.1 completed!\n$(tput sgr 0)"

exit 0
