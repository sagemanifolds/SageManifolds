#!/bin/bash

#cd $SAGE_ROOT

echo -e "\n$(tput setaf 4)Downloading the package SageManifolds v0.5...$(tput sgr 0)"
# downloading with either wget or curl 
wget -N http://sagemanifolds.obspm.fr/spkg/manifolds-0.5.tar.gz || curl -O http://sagemanifolds.obspm.fr/spkg/manifolds-0.5.tar.gz

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

#
# This part modifies the default output of derivatives
#
echo -e "\n$(tput setaf 4)Downloading modified pynac functions (pynac_mod.pxi)...$(tput sgr 0)"
# downloading with either wget or curl 
wget -N http://sagemanifolds.obspm.fr/spkg/pynac_mod.pxi || curl -O http://sagemanifolds.obspm.fr/spkg/pynac_mod.pxi

if [ ! -f "pynac_mod.pxi" ]; then
	echo  -e "\n$(tput setaf 1)Download the pynac_mod.pxi manually and run again$(tput sgr 0)"
	exit 0 
else 
	mv pynac_mod.pxi src/sage/symbolic/
fi

# do a backup of pynac.pyx 'just in case' 
if [ ! -f "src/sage/symbolic/pynac.pyx_orig" ]; then
	cp src/sage/symbolic/pynac.pyx src/sage/symbolic/pynac.pyx_orig
fi 

func=('def py_print_function_pystring' 
		'def py_latex_function_pystring'
		'cdef public stdstring\* py_print_fderivative' 
		'cdef public stdstring\* py_latex_fderivative')

echo -e "$(tput setaf 4)In src/sage/symbolic/pynac.pyx:$(tput sgr 0)"
let i=0
while (( ${#func[@]} > i )); do
    echo -e "$(tput setaf 4)Commenting out ${func[i]}...$(tput sgr 0)"
	sed -i '/^'"${func[i]}"'/,/^def\|^cdef/{/^def\|^cdef/!s/^/#/g};/^'"${func[i]}"'/s/^/#/' src/sage/symbolic/pynac.pyx
	((i++))
done

if [[ -z $(grep "pynac_mod.pxi" src/sage/symbolic/pynac.pyx) ]]; then
	echo -e "$(tput setaf 4)Including the pynac_mod.pxi file...$(tput sgr 0)"
	echo 'include "sage/symbolic/pynac_mod.pxi"' >> src/sage/symbolic/pynac.pyx
fi 

#
# Re-building sage  
#
echo -e "\n$(tput bold)$(tput setaf 4)Running ./sage -b to re-build sage.$(tput sgr 0)"
./sage -b 

echo -e "\n$(tput bold)$(tput setaf 4)Building the html documentation (may take some time).$(tput sgr 0)"
./sage -docbuild reference inventory
./sage -docbuild reference html
./sage -docbuild tensors_free_module html
./sage -docbuild sagemanifolds html

echo -e "\n$(tput bold)$(tput setaf 4)Installation of SageManifolds completed!$(tput sgr 0)"

exit 0

