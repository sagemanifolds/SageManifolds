#!/bin/bash

echo -e "\n$(tput setaf 4)Downloading modified pynac functions...$(tput sgr 0)"
# downloading with either wget or curl 
wget -N http://sagemanifolds.obspm.fr/spkg/pynac_mod.pxi || curl -O http://sagemanifolds.obspm.fr/spkg/pyna_mod.pxi

if [ ! -f "pynac_mod.pxi" ]; then
	echo  -e "\n$(tput setaf 1)Download the file manually and run again$(tput sgr 0)"
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

let i=0
while (( ${#func[@]} > i )); do
    printf "Commenting out ${func[i]}...\n"
	sed -i '/^'"${func[i]}"'/,/^def\|^cdef/{/^def\|^cdef/!s/^/#/g};/^'"${func[i]}"'/s/^/#/' src/sage/symbolic/pynac.pyx
	((i++))
done

echo 'include "sage/symbolic/pynac_mod.pxi"' >> src/sage/symbolic/pynac.pyx

exit 0
