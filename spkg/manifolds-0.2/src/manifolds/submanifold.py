r"""
Submanifolds

The class :class:`Submanifold` implements submanifolds of differentiable 
manifolds over `\RR`. By 'submanifold' it is meant embedded manifold. 

A subclass of :class:`Submanifold` is naturally :class:`MCurve` (curves as 
1-dimensional submanifolds). 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES:
    
    The sphere `S^2` as a submanifold of `\RR^3`::
    
        sage: m = Manifold(3, 'R3', r'\mathcal{M}')
        sage: c_cart = Chart(m, 'x y z', 'cart') # Cartesian coordinates set on R^3 
        sage: (u,v) = var('u v')  # spherical coordinates on the sphere
        sage: s = Submanifold(m, 2, 'u:[0,pi] v:[0,2*pi)', 'spher', [sin(u)*cos(v), sin(u)*sin(v), cos(u)], 'S2', r'\mathcal{S}')
        sage: s
        2-dimensional manifold 'S2'
        sage: s.embedding
        differentiable mapping from 2-dimensional manifold 'S2' to 3-dimensional manifold 'R3'
        sage: s.embedding.coord_expression
        {('spher', 'cart'): functions (cos(v)*sin(u), sin(u)*sin(v), cos(u)) on the chart 'spher' (u, v)}
        
    A helix as a submanifold of `\RR^3`::
    
        sage: t = var('t') # parameter along the helix
        sage: h = Submanifold(m, 1, 't', 't', [cos(t), sin(t), t], 'helix', 'H')
        sage: h
        1-dimensional manifold 'helix'
        
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
from manifold import Manifold
from chart import Chart
from diffmapping import DiffMapping

class Submanifold(Manifold):
    r"""
    Base class for submanifolds, i.e. manifolds embedded in a differentiable 
    manifold.
    
    INPUT:
    
    - ``ambient_manifold`` -- the ambient manifold
    - ``n`` -- dimension of the submanifold
    - ``coordinates`` -- single string defining the coordinate symbols and 
      ranges of a chart defined on the submanifold: see the documentation of
      :class:`Chart` for the syntax
    - ``chart_name`` -- string containing the name given to the chart 
      corresponding to the above coordinates (should be rather short)
    - ``embedding_functions`` -- the coordinate expression of the embedding in 
      the chart defined above, the arrival chart being ``ambient_chart``: 
      list of symbolic expressions (or strings if the involved symbols have not
      been previously defined), each item representing a coordinate 
      expression
    - ``name`` -- (default: None) name given to the submanifold 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the submanifold
    - ``ambient_chart`` -- (default: None) the chart of the ambient manifold 
      for defining the embedding; if none is provided the ambient manifold's
      default chart will be used
    - ``start_index`` -- (default: 0) lower bound of the range of indices on the
      submanifold

    EXAMPLES:
        
        The sphere `S^2` as a submanifold of `\RR^3`::
    
            sage: m = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart = Chart(m, 'x y z', 'cart') # Cartesian coordinates set on R^3 
            sage: (u,v) = var('u v')  # spherical coordinates on the sphere
            sage: s = Submanifold(m, 2, 'u:[0,pi] v:[0,2*pi)', 'spher', [sin(u)*cos(v), sin(u)*sin(v), cos(u)], 'S2', r'\mathcal{S}')
            sage: s
            2-dimensional manifold 'S2'
            sage: s.embedding
            differentiable mapping from 2-dimensional manifold 'S2' to 3-dimensional manifold 'R3'
            sage: s.embedding.coord_expression
            {('spher', 'cart'): functions (cos(v)*sin(u), sin(u)*sin(v), cos(u)) on the chart 'spher' (u, v)}
        
        The coordinate expressions of the embedding can be passed as strings if
         the corresponding symbols have not been previously defined::
          
            sage: h2 = Submanifold(m, 2, 'r:[0,+oo) ph:[0,2*pi)', 'spher', ["sinh(r)*cos(ph)", "sinh(r)*sin(ph)", "cosh(r)"], 'H2')
            sage: h2.def_chart
            chart 'spher' (r, ph)
            

 
    """
    def __init__(self, ambient_manifold, n, coordinates, chart_name, 
                 embedding_functions, name=None, latex_name=None, 
                 ambient_chart=None, start_index=0):
        from sage.symbolic.ring import SR
        Manifold.__init__(self, n, name, latex_name, start_index)
        if not isinstance(ambient_manifold, Manifold):
            raise TypeError("The argument ambient_manifold must be a manifold.")
        self.ambient_manifold = ambient_manifold
        Chart(self, coordinates, chart_name)
        if ambient_chart is None:
            ambient_chart = ambient_manifold.def_chart.name
        n_amb = ambient_manifold.dim
        if len(embedding_functions) != n_amb:
            raise ValueError(str(n_amb) + 
                             " coordinate functions must be provided.")
        embedding_expressions = [SR(embedding_functions[i]) for i in
                                 range(n_amb)]
        self.embedding = DiffMapping(self, ambient_manifold, embedding_expressions, 
                                     chart_name, ambient_chart)


    def plot(self, coord_ranges, chartname = None, **kwds):
        r"""
        Plot of a submanifold embedded in `\RR^2` or `\RR^3`

        INPUT:

         - ``coord_ranges`` -- list of pairs (u_min, u_max) for each coordinate
           u on the submanifold
         - ``chartname`` -- (default: None) name of the chart in which the 
           above coordinates are defined; if none is provided, the submanifold
           default chart is assumed. 
         - ``**kwds`` -- (default: None) keywords passed to Sage graphic 
           routines
           
        OUTPUT:
        
          - Graphics3d object (ambient manifold = `\RR^3`) or Graphics object
            (ambient manifold = `\RR^2`)

        EXAMPLES:

            Plot of a torus embedded in `\RR^3`::
                
                sage: forget() # required for the doctest (to forget previous assumptions)
                sage: m = Manifold(3, 'R3')
                sage: c_cart  = Chart(m, 'x y z', 'cart')
                sage: (u,v) = var('u v')
                sage: t = Submanifold(m, 2, 'u:[0,2*pi) v:[0,2*pi)', "canonical", [(2+cos(u))*cos(v),(2+cos(u))*sin(v),sin(u)], "torus")
                sage: t.plot([[0,2*pi], [0,2*pi]], aspect_ratio=1)


            Plot of a helix embedded in `\RR^3`::
                
                sage: h = Submanifold(m, 1, 'u', 'u', [cos(u), sin(u), u], "helix")
                sage: h.plot([0,20])
                
               
            Plot of an Archimedean spiral embedded in `\RR^2`::
            
                sage: R2 = Manifold(2, 'R2')
                sage: c_cart = Chart(R2, 'x y', 'cart')
                sage: t = var('t')
                sage: s = Submanifold(R2, 1, 't', 't', [t*cos(t), t*sin(t)], "spiral")
                sage: s.plot([0,40])

        """
        from sage.plot.plot import parametric_plot
        if chartname is None:
            chartname = self.def_chart.name
        amb = self.ambient_manifold.name
        if amb == 'R3' or amb == 'R2':
            if 'cart' not in self.ambient_manifold.atlas:
                raise ValueError("For drawing, " + amb + 
                                 " Cartesian coordinates must be defined.")
            coord_functions = \
                self.embedding.coord_expression[(chartname, 'cart')].functions
            chart = self.atlas[chartname]
            if self.dim == 1:
                urange = (chart.xx[0], coord_ranges[0], coord_ranges[1])
                graph = parametric_plot(coord_functions, urange, **kwds)
            elif self.dim == 2: 
                urange = (chart.xx[0], coord_ranges[0][0], coord_ranges[0][1])
                vrange = (chart.xx[1], coord_ranges[1][0], coord_ranges[1][1])
                graph = parametric_plot(coord_functions, urange, vrange, **kwds)
            else: 
                raise ValueError("The dimension must be at most 2 " + 
                                 "for plotting.")
        else:
            raise NotImplementedError("Plotting is implemented only for " + 
                                      "submanifolds of R^2 or R^3.")
        return graph


#******************************************************************************

class MCurve(Submanifold):
    r"""
    Class for curves, i.e. 1-dimensional submanifolds embedded in a 
    differentiable manifold.
    
    
    .. NOTE::

        Ideally the name should be Curve instead of MCurve, but this would
        clash with the function 
        :func:`sage.schemes.plane_curves.constructor.Curve`.
    
    
    INPUT:
    
    - ``ambient_manifold`` -- the manifold in which the curve is defined
    - ``paramater`` -- string defining the parameter along the curve; it has 
       at most three fields, separated by ':': 
        
        1. The parameter symbol (a letter or a few letters)
        2. (optional) The interval defining the parameter range: if not
           provided, the parameter is assumed to span all `\RR`; otherwise 
           the interval must be provided in the form [a,b], (a,b] (or 
           equivalently ]a,b]), (a,b) (or equivalently ]a,b[).
           The bounds a and b can be +/-Infinity, Inf, infinity, inf or oo.
           Note that the interval declaration must not contain any space 
           character.
        3. (optional) The LaTeX spelling of the parameter; if not provided the
           coordinate symbol given in the first field will be used.

      The order of the fields 2 and 3 does not matter and each of them can be
      omitted.
      If it contains any LaTeX expression, the string ``parameter`` must be
      declared with the prefix 'r' (for "raw") to allow for a proper treatment 
      of the backslash character. 

    - ``param_name`` -- string for the name given to the curve parametrization
      (actually name of the chart on the curve manifold)
    - ``embedding_functions`` -- the coordinate expression of the embedding in 
      terms of the parameter defined above, the arrival chart being 
      ``ambient_chart`` : list of symbolic expressions (or strings if the 
      involved symbols have not been previously defined), each item 
      representing a coordinate expression
    - ``name`` -- (default: None) name given to the curve
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the curve
    - ``ambient_chart`` -- (default: None) the chart of the ambient manifold 
      for defining the curve embedding; if none is provided the ambient 
      manifold's default chart will be used

    EXAMPLES:
        
        A helix in `\RR^3`::
        
            sage: m = Manifold(3, 'R3')
            sage: c_cart  = Chart(m, 'x y z', 'cart')
            sage: t = var('t')
            sage: h = MCurve('t','t',m, [cos(t), sin(t), t], "helix")
            sage: h
            curve 'helix' on 3-dimensional manifold 'R3'
            
        A curve is a manifold::
        
            sage: isinstance(h, Manifold)
            True
            sage: h.atlas
            {'t': chart 't' (t,)}
            sage: p = Point(h, (pi,)) ; p
            point on curve 'helix' on 3-dimensional manifold 'R3'
            
        A curve is a submanifold::
        
            sage: p1 = h.embedding(p) ; p1
            point on 3-dimensional manifold 'R3'
            sage: p1.coord()
            (-1, 0, pi)
            sage: h.plot([0,20])
 
    """
    def __init__(self, parameter, param_name, ambient_manifold,
                 embedding_functions, name=None, latex_name=None, 
                 ambient_chart=None):
        Submanifold.__init__(self, ambient_manifold, 1, parameter, 
                             param_name, embedding_functions, name, 
                             latex_name, ambient_chart)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "curve"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on " + str(self.ambient_manifold)
        return description

