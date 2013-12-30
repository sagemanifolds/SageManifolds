r"""
Vector fields

The class :class:`VectorField` implements vector fields on differentiable 
manifolds over `\RR`. 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

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

from tensorfield import TensorField
from scalarfield import ScalarField

class VectorField(TensorField) :
    r"""
    Class for vector fields on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the vector field is defined
      (must be an instance of class :class:`Domain`)
    - ``name`` -- (default: None) name given to the vector field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the vector field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = m.chart('x y z', 'xyz-coord')
        sage: v = VectorField(m, 'V') ; v
        vector field 'V' on the 3-dimensional manifold 'M'
        sage: latex(v)
        V

    A vector field is a tensor field of rank 1 and of type (1,0)::
    
        sage: v.rank
        1
        sage: v.tensor_type
        (1, 0)

    Components of a vector field with respect to a given frame::
    
        sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
        sage: v[0], v[1], v[2] = (1, 4, 9)  # components on M's default frame (e)
        sage: v.comp() 
        1-index components w.r.t. the vector frame 'e'
    
    The totality of the components are accessed via the operator [:]::
    
        sage: v[:] = (1, 4, 9) # equivalent to v[0], v[1], v[2] = (1, 4, 9)
        sage: v[:]
        [1, 4, 9]
        
    The components are also read on the expansion on the frame 'e', as provided
    by the method :meth:`view`::
    
        sage: v.view()   # displays the expansion on the manifold's default frame (e)
        V = e_0 + 4 e_1 + 9 e_2
    
    A subset of the components can be accessed by means of Python's slice 
    notation::
        
        sage: v[1:] = (-2, -3)
        sage: v[:]
        [1, -2, -3]
        sage: v[:2]
        [1, -2]
        
    The components are instances of the class :class:`Components`::
    
        sage: type(v.comp())  
        <class 'sage.geometry.manifolds.component.Components'>

    Components in another frame::
    
        sage: f = VectorFrame(m, 'f')
        sage: for i in range(3):
        ...       v.set_comp('f')[i] = (i+1)**3
        ...
        sage: v.comp('f')[2]
        27
        sage: v.view('f')
        V = f_0 + 8 f_1 + 27 f_2

    The range of the indices depends on the convention set for the manifold::
        
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = m.chart('x y z', 'xyz-coord')
        sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
        sage: v = VectorField(m, 'V')
        sage: (v[1], v[2], v[3]) = (1, 4, 9)
        sage: v[0]
        Traceback (most recent call last):
        ...
        IndexError: Index out of range: 0 not in [1,3]

    A vector field acts on scalar fields (derivation along the vector field)::
    
        sage: m = Manifold(2, 'M')            
        sage: c_cart.<x,y> = m.chart('x y', 'cart')
        sage: f = ScalarField(m, x*y^2, name='f')  
        sage: v = VectorField(m, 'v')         
        sage: v[:] = (-y, x)
        sage: v.view()
        v = -y d/dx + x d/dy
        sage: v(f)
        scalar field 'v(f)' on the 2-dimensional manifold 'M'
        sage: v(f).expr()
        2*x^2*y - y^3
        sage: latex(v(f))
        v\left(f\right)

    """
    def __init__(self, domain, name=None, latex_name=None) :
        TensorField.__init__(self, domain, 1, 0, name, latex_name)
        self._init_dependencies()
        
    def _repr_(self) :
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "vector field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description


    def _new_instance(self):
        r"""
        Create a :class:`VectorField` instance. 
        
        """
        return VectorField(self.domain)

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        TensorField._del_derived(self)
        self._del_dependencies()
        
    def _init_dependencies(self):
        r"""
        Initialize list of quantities that depend on ``self``
        """
        self._lie_der_along_self = {}

    def _del_dependencies(self):
        r"""
        Clear list of quantities that depend on ``self``
        """
        if self._lie_der_along_self != {}:
            for idtens, tens in self._lie_der_along_self.items():
                del tens._lie_derivatives[id(self)]
            self._lie_der_along_self.clear()

        
    def __call__(self, scalar):
        r"""
        Action on a scalar field.
            
        INPUT:
            
        - ``scalar`` -- scalar field `f`
            
        OUTPUT:
            
        - scalar field formed by the derivative of `f` along the vector 
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`
          
        EXAMPLES:
        
        Action of a vector field on a scalar field on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')            
            sage: c_cart.<x,y> = m.chart('x y', 'cart')
            sage: f = ScalarField(m, x*y^2)  
            sage: v = VectorField(m, 'v')         
            sage: v[:] = (-y, x)
            sage: v(f)
            scalar field on the 2-dimensional manifold 'M'
            sage: v(f).expr()
            2*x^2*y - y^3
          
        """
        from diffform import OneForm
        from scalarfield import ZeroScalarField
        if isinstance(scalar, OneForm):
            # This is actually the action of the vector field on a 1-form, 
            # as a tensor field of type (1,0):
            return scalar(self)
        if not isinstance(scalar, ScalarField):
            raise TypeError("The argument must be a scalar field")
        if scalar.manifold != self.manifold:
            raise ValueError("The scalar field and the vector field must" + 
                                 " be defined on the same manifold.")
        if isinstance(scalar, ZeroScalarField):
            return scalar
        # search for a commont chart: 
        chart_name = None
        def_chart_name = self.domain.def_chart.name
        if def_chart_name in scalar.express:
            if def_chart_name + '_b' in self.components:
                chart_name = def_chart_name
        else:
            for kchart in scalar.express:
                if kchart + '_b' in self.components: 
                    chart_name = kchart
                    break
        if chart_name is None:
            raise ValueError("No common chart found.")
        chart = self.domain.atlas[chart_name]
        v = self.comp(chart_name + "_b")
        f = scalar.function_chart(chart_name) 
        si = self.manifold.sindex
        res = 0 
        for i in range(self.manifold.dim):
            res += v[i+si, chart_name] * f.diff(chart.xx[i])
        # Name of the output:
        res_name = None
        if self.name is not None and scalar.name is not None:
            res_name = self.name + "(" + scalar.name + ")"
        # LaTeX symbol for the output:
        res_latex = None
        if self.latex_name is not None and scalar.latex_name is not None:
            res_latex = self.latex_name + r"\left(" + scalar.latex_name + \
                        r"\right)"
        return res.scalar_field(name=res_name, latex_name=res_latex)

