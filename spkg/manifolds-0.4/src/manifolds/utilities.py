r"""
SageManifolds utilities. 

This module defines helper functions that are not class methods.


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version
- Joris Vankerschaver (2010): for the function is_atomic()

"""

#******************************************************************************
#       Copyright (C) 2013 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject

def is_atomic(expression):
    r"""
    Helper function to check whether some LaTeX expression is atomic.
    
    Adapted from function :meth:`DifferentialFormFormatter._is_atomic` written 
    by Joris Vankerschaver (2010)
    
    INPUT:
    
    - ``expression`` -- string representing the expression (e.g. LaTeX string)
    
    OUTPUT:
    
    - True if additive operations are enclosed in parentheses, false otherwise.

    EXAMPLES::
    
        sage: from sage.geometry.manifolds.utilities import is_atomic
        sage: is_atomic("2*x")
        True
        sage: is_atomic("2+x")
        False
        sage: is_atomic("(2+x)")
        True
        sage: m = Manifold(2, 'M')
        sage: c_xy.<x,y> = m.chart('x y')
        sage: f = FunctionChart(c_xy, x^2+3*y+1) ; f
        x^2 + 3*y + 1
        sage: is_atomic(latex(f))
        False
        sage: g = FunctionChart(c_xy, 3*x*y) ; g
        3*x*y
        sage: is_atomic(latex(g))
        True

    """
    if not isinstance(expression, basestring):
        raise TypeError("The argument must be a string.")
    level = 0
    for n, c in enumerate(expression):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        if c == '+' or c == '-':
            if level == 0 and n > 0:
                return False
    return True


def is_atomic_wedge_txt(expression):
    r"""
    Helper function to check whether some text-formatted expression is atomic 
    in terms of wedge products. 
    
    Adapted from function :meth:`DifferentialFormFormatter._is_atomic` written 
    by Joris Vankerschaver (2010)
    
    INPUT:
    
    - ``expression`` -- string representing the text-formatted expression
    
    OUTPUT:
    
    - True if wedge products are enclosed in parentheses, false otherwise.

    EXAMPLES::
    
        sage: from sage.geometry.manifolds.utilities import is_atomic_wedge_txt
        sage: is_atomic_wedge_txt("a")
        True
        sage: is_atomic_wedge_txt(r"a/\b")
        False
        sage: is_atomic_wedge_txt(r"(a/\b)")
        True
        sage: is_atomic_wedge_txt(r"(a/\b)/\c")
        False
        sage: is_atomic_wedge_txt(r"(a/\b/\c)")
        True

    """
    if not isinstance(expression, basestring):
        raise TypeError("The argument must be a string.")
    level = 0
    for n, c in enumerate(expression):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        if c == '/' and expression[n+1:n+2] == '\\':
            if level == 0 and n > 0:
                return False
    return True


def is_atomic_wedge_latex(expression):
    r"""
    Helper function to check whether LaTeX-formatted expression is atomic in 
    terms of wedge products.
    
    Adapted from function :meth:`DifferentialFormFormatter._is_atomic` written 
    by Joris Vankerschaver (2010)
    
    INPUT:
    
    - ``expression`` -- string representing the LaTeX expression
    
    OUTPUT:
    
    - True if wedge products are enclosed in parentheses, false otherwise.

    EXAMPLES::
    
        sage: from sage.geometry.manifolds.utilities import is_atomic_wedge_latex
        sage: is_atomic_wedge_latex(r"a")                    
        True
        sage: is_atomic_wedge_latex(r"a\wedge b")
        False
        sage: is_atomic_wedge_latex(r"(a\wedge b)")
        True
        sage: is_atomic_wedge_latex(r"(a\wedge b)\wedge c")
        False
        sage: is_atomic_wedge_latex(r"((a\wedge b)\wedge c)")
        True
        sage: is_atomic_wedge_latex(r"(a\wedge b\wedge c)")  
        True
        sage: is_atomic_wedge_latex(r"\omega\wedge\theta")      
        False
        sage: is_atomic_wedge_latex(r"(\omega\wedge\theta)")
        True
        sage: is_atomic_wedge_latex(r"\omega\wedge(\theta+a)")
        False

    """
    if not isinstance(expression, basestring):
        raise TypeError("The argument must be a string.")
    level = 0
    for n, c in enumerate(expression):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        if c == '\\' and expression[n+1:n+6] == 'wedge':
            if level == 0 and n > 0:
                return False
    return True


def format_mul_txt(name1, operator, name2):
    r"""
    Helper function for text-formatted names of results of multiplication or 
    tensor product. 
    
    """
    if name1 is None or name2 is None:
        return None
    if not is_atomic(name1) or not is_atomic_wedge_txt(name1):
        name1 = '(' + name1 + ')'
    if not is_atomic(name2) or not is_atomic_wedge_txt(name2):
        name2 = '(' + name2 + ')'
    return name1 + operator + name2 


def format_mul_latex(name1, operator, name2):
    r"""
    Helper function for LaTeX names of results of multiplication or tensor 
    product. 
    
    """
    if name1 is None or name2 is None:
        return None
    if not is_atomic(name1) or not is_atomic_wedge_latex(name1):
        name1 = r'\left(' + name1 + r'\right)'
    if not is_atomic(name2) or not is_atomic_wedge_latex(name2):
        name2 = r'\left(' + name2 + r'\right)'
    return name1 + operator + name2 


def format_unop_txt(operator, name):
    r"""
    Helper function for text-formatted names of results of unary operator.
    
    """
    if name is None:
        return None
    if not is_atomic(name) or not is_atomic_wedge_txt(name):
    #!# is_atomic_otimes_txt should be added
        name = '(' + name + ')'
    return operator + name


def format_unop_latex(operator, name):
    r"""
    Helper function for LaTeX names of results of unary operator.
    
    """
    if name is None:
        return None
    if not is_atomic(name) or not is_atomic_wedge_latex(name):
    #!# is_atomic_otimes_latex should be added
        name = r'\left(' + name + r'\right)'
    return operator + name

class FormattedExpansion(SageObject):
    r""" 
    Helper class for displaying tensor expansions.
    """
    def  __init__(self, tensor):
        self.tensor = tensor
        self.txt = None
        self.latex = None
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return self.txt
        
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return self.latex


#***********************************************************

def simple_determinant(aa):
    r"""
    Compute the determinant of a square matrix.
    
    This function is a workaround to bypass a bug in Sage det method.
    """
    from sage.matrix.constructor import matrix
    n = aa.nrows()
    if n == 1:
        return aa[0,0]
    res = 0
    sign = True
    for i in range(n):
        b = []
        for k in range(i):
            r = []
            for l in range(1,n):
               r.append(aa[k,l])
            b.append(r)
        for k in range(i+1,n):
            r = []
            for l in range(1,n):
               r.append(aa[k,l])
            b.append(r)
        bb = matrix(b)
        if sign:
            res += aa[i,0] * simple_determinant(bb)
        else:
            res -= aa[i,0] * simple_determinant(bb)
        sign = not sign
    return res



def simplify_sqrt_real(expr):
    r"""
    Simplify sqrt in symbolic expressions in the real domain.
    
    EXAMPLES:
    
    Simplifications of basic expressions::
    
        sage: assume(x<0)      
        sage: simplify_sqrt_real( sqrt(x^2) )
        -x
        sage: simplify_sqrt_real( sqrt(x^2-2*x+1) )
        -x + 1
        sage: simplify_sqrt_real( sqrt(x^2) + sqrt(x^2-2*x+1) )
        -2*x + 1

    """
    from sage.symbolic.ring import SR
    from sage.calculus.calculus import maxima
    # 1/ Search for the sqrt's in expr
    sexpr = str(expr)
    if 'sqrt(' not in sexpr:  # no sqrt to simplify
        return expr
    pos_sqrts = []   # positions of the sqrt's in sexpr
    the_sqrts = []   # the sqrt sub-expressions in sexpr
    for pos in range(len(sexpr)):
        if sexpr[pos:pos+5] == 'sqrt(':
            pos_sqrts.append(pos)
            parenth = 1
            scan = pos+5
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1 
            the_sqrts.append( sexpr[pos:scan] )
    # 2/ Simplifications of the sqrt's
    new_expr = ""    # will contain the result
    pos0 = 0
    for i, pos in enumerate(pos_sqrts):
        # radcan is called on each sqrt:
        x = SR(the_sqrts[i])
        simpl = SR(x._maxima_().radcan())
        # the absolute value of radcan's output is taken, the call to simplify() 
        # taking into account possible assumptions regarding the sign of simpl:
        new_expr += sexpr[pos0:pos] + '(' + str(abs(simpl).simplify()) + ')'
        pos0 = pos + len(the_sqrts[i])
    new_expr += sexpr[pos0:]
    return SR(new_expr)

def simplify_abs_trig(expr):
    r"""
    Simplify abs(sin(...)) in symbolic expressions
    """
    from sage.symbolic.ring import SR
    from sage.symbolic.constants import pi
    sexpr = str(expr)
    if 'abs(sin(' not in sexpr:  # nothing to simplify
        return expr
    tp = []
    val = []
    for pos in range(len(sexpr)):
        if sexpr[pos:pos+8] == 'abs(sin(':
            # finding the end of abs argument:
            scan = pos+4 # start of abs
            parenth = 1
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1 
            pos_abs_end = scan
            # finding the end of sin argument:
            scan = pos+8 # start of sin
            parenth = 1
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1 
            pos_sin_end = scan
            # if the abs contains only the sinus, the simplification can be tried:
            if pos_sin_end == pos_abs_end-1: 
                tp.append(pos)
                val.append( sexpr[pos:pos_abs_end] )
    simp = []
    for v in val:
        # argument of the sinus:
        sx = v[8:-2]
        x = SR(sx)
        if x>=0 and x<=pi:
            simp.append('sin(' + sx + ')')
        elif x>=-pi and x<=0:
            simp.append('(-sin(' + sx + '))')
        else:
            simp.append(v)  # no simplification is applicable
    nexpr = ""
    pos0 = 0
    for i, pos in enumerate(tp):
        nexpr += sexpr[pos0:pos] + simp[i] 
        pos0 = pos + len(val[i])
    nexpr += sexpr[pos0:]
    return SR(nexpr)
    


def simplify_chain(expr):
    r"""
    Perform a chain of simplications to a symbolic expression.
    
    """
    expr = expr.simplify_factorial()
    expr = expr.simplify_trig()
    expr = expr.simplify_rational()
    expr = simplify_sqrt_real(expr)
    expr = simplify_abs_trig(expr)
    expr = expr.simplify_radical()
    expr = expr.simplify_log('one')
    expr = expr.simplify_rational()
    expr = expr.simplify_trig()
    return expr

    
    
