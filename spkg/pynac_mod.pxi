def py_print_function_pystring(id, args, fname_paren=False):
    """
    Return a string with the representation of the symbolic function specified
    by the given id applied to args.

    INPUT:

        id --   serial number of the corresponding symbolic function
        params -- Set of parameter numbers with respect to which to take
                    the derivative.
        args -- arguments of the function.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_print_function_pystring, get_ginac_serial
        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: var('x,y,z')
        (x, y, z)
        sage: foo = function('foo', nargs=2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_function_pystring(i, (x,y))
        'foo(x, y)'
        sage: py_print_function_pystring(i, (x,y), True)
        '(foo)(x, y)'
        sage: def my_print(self, *args): return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('foo', nargs=2, print_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_function_pystring(i, (x,y))
        'my args are: x, y'
    """
    cdef Function func = get_sfunction_from_serial(id)
    # This function is called from two places, from function::print in Pynac
    # and from py_print_fderivative. function::print checks if the serial
    # belongs to a function defined at the C++ level. There are no C++ level
    # functions that return an instance of fderivative when derivated. Hence,
    # func will never be None.
    assert(func is not None)

    # if function has a custom print function call it
    if hasattr(func,'_print_'):
        res = func._print_(*args)
        # make sure the output is a string
        if res is None:
            return ""
        if not isinstance(res, str):
            return str(res)
        return res

    # otherwise use default output
    if fname_paren:
        olist = ['(', func._name, ')']
    else:
        olist = [func._name]

    # default: print the arguments
    olistnp = ''.join(olist)
    olist.extend(['(', ', '.join(map(repr, args)), ')'])

    try:
        omit_function_args_choice
    except NameError:
        pass
    else:
        if omit_function_args_choice:
            olist = olistnp 

    return ''.join(olist)

def py_latex_function_pystring(id, args, fname_paren=False):
    """
    Return a string with the latex representation of the symbolic function
    specified by the given id applied to args.

    See documentation of py_print_function_pystring for more information.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_function_pystring, get_ginac_serial
        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: var('x,y,z')
        (x, y, z)
        sage: foo = function('foo', nargs=2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        '{\\rm foo}\\left(x, y^{z}\\right)'
        sage: py_latex_function_pystring(i, (x,y^z), True)
        '\\left({\\rm foo}\\right)\\left(x, y^{z}\\right)'
        sage: py_latex_function_pystring(i, (int(0),x))
        '{\\rm foo}\\left(0, x\\right)'

    Test latex_name::

        sage: foo = function('foo', nargs=2, latex_name=r'\mathrm{bar}')
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        '\\mathrm{bar}\\left(x, y^{z}\\right)'

    Test custom func::

        sage: def my_print(self, *args): return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('foo', nargs=2, print_latex_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+100):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        'my args are: x, y^z'


    """
    cdef Function func = get_sfunction_from_serial(id)
    # This function is called from two places, from function::print in Pynac
    # and from py_latex_fderivative. function::print checks if the serial
    # belongs to a function defined at the C++ level. There are no C++ level
    # functions that return an instance of fderivative when derivated. Hence,
    # func will never be None.
    assert(func is not None)

    # if function has a custom print method call it
    if hasattr(func, '_print_latex_'):
        res = func._print_latex_(*args)
        # make sure the output is a string
        if res is None:
            return ""
        if not isinstance(res, str):
            return str(res)
        return res

    # otherwise, use the latex name if defined
    if func._latex_name:
        name = func._latex_name
    else:
        # if latex_name is not defined, then call
        # latex_variable_name with "is_fname=True" flag
        from sage.misc.latex import latex_variable_name
        name = latex_variable_name(func._name, is_fname=True)
    if fname_paren:
        olist = [r'\left(', name, r'\right)']
    else:
        olist = [name]

    # default: print the arguments
    olistnp = ''.join(olist)
    from sage.misc.latex import latex
    olist.extend([r'\left(', ', '.join([latex(x) for x in args]), r'\right)'] )

    try:
        omit_function_args_choice
    except NameError:
        pass
    else:
        if omit_function_args_choice:
            olist = olistnp

    return ''.join(olist)

cdef public stdstring* py_print_fderivative(unsigned id, object params,
        object args) except +:
    """
    Return a string with the representation of the derivative of the symbolic
    function specified by the given id, lists of params and args.

    INPUT:

        id --   serial number of the corresponding symbolic function
        params -- Set of parameter numbers with respect to which to take
                    the derivative.
        args -- arguments of the function.


    """
    ostr = ''.join(['D[', ', '.join([repr(int(x)) for x in params]), ']'])
    bra  = True

    try: 
        textbook_style_deriv_choice
    except NameError:
        pass
    else:
        if textbook_style_deriv_choice:
            if(len(params)>1):
                op = ''.join(['D^',str(len(params)),'/D'])
            else: 
                op = 'D/D' 
            ostr = ''.join([op, 'D'.join([repr(args[int(x)]) for x in params]), ' '])
            bra = False

    fstr = py_print_function_pystring(id, args, bra)
    py_res = ostr + fstr
    return string_from_pystr(py_res)

def textbook_style_deriv(c=False):
        global textbook_style_deriv_choice
        textbook_style_deriv_choice = c

def omit_function_args(c=False):
        global omit_function_args_choice
        omit_function_args_choice = c

cdef public stdstring* py_latex_fderivative(unsigned id, object params,
        object args) except +:
    """
    Return a string with the latex representation of the derivative of the
    symbolic function specified by the given id, lists of params and args.

    See documentation of py_print_fderivative for more information.

    """
    ostr = ''.join(['D[', ', '.join([repr(int(x)) for x in params]), ']'])
    bra  = True

    try:
        textbook_style_deriv_choice
    except NameError:
        pass
    else:
        if textbook_style_deriv_choice:
            from sage.misc.latex import latex 
            if(len(params)>1):
                op = ''.join(['\\frac{\partial^',str(len(params)),'}{\partial '])
            else: 
                op = '\\frac{\partial}{\partial '

            ostr = ''.join([op, '\partial '.join([''.join([latex(args[int(x)]), '^{', str(params.count(int(x)) if  params.count(int(x)) > 1 else ""), '}']) for x in list(set(params))]), '}'])


            bra = False

    fstr = py_latex_function_pystring(id, args, bra)
    py_res = ostr + fstr
    return string_from_pystr(py_res)
