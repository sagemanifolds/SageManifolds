r"""
Differential forms

The class :class:`DiffForm` implements differential forms on differentiable 
manifolds over `\RR`. 

It is a subclass of :class:`TensorField`, differential forms being a special 
type of tensor fields. 

Subclasses of :class:`DiffForm` are

* :class:`ScalarField` for differential forms of degree 0 (i.e. scalar fields)
* :class:`OneForm` for differential forms of degree 1 (i.e. 1-forms). 

.. NOTE::

    A difference with the preceding Sage class :class:`DifferentialForm` 
    is that the present class lies at the tensor field level. Accordingly, an
    instance of :class:`DiffForm` can have various sets of components, each in
    a different coordinate system or coframe, while the class 
    :class:`DifferentialForm` considers differential forms at the component 
    level in a fixed chart. In this respect, the class 
    :class:`DifferentialForm` is closer to the class 
    :class:`CompFullyAntiSym` than to :class:`DiffForm`

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version
- Joris Vankerschaver (2010): developed a previous class, 
  :class:`DifferentialForm` (cf. the above note), which inspired the storage of 
  the non-zero components as a dictionary whose keys are the indices.

"""

#******************************************************************************
#       Copyright (C) 2013 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2010 Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from tensorfield import TensorField
from component import CompFullyAntiSym


class DiffForm(TensorField):
    r"""
    Class for differential forms on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the differential form is 
      defined (must be an instance of class :class:`Domain`)
    - ``p`` -- the degree of the differential form (i.e. its tensor rank)
    - ``name`` -- (default: None) name given to the differential form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the differential 
      form; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    A 2-form on a 4-dimensional manifold::
    
        sage: m = Manifold(4, 'M')
        sage: c_txyz = m.chart('t x y z')
        sage: a = DiffForm(m, 2, 'a') ; a
        2-form 'a' on the 4-dimensional manifold 'M'
        
    A differential form is a tensor field of purely covariant type::
    
        sage: isinstance(a, TensorField)
        True
        sage: a.tensor_type  
        (0, 2)

    It is antisymmetric, its components being instances of the class 
    :class:`CompFullyAntiSym`::
    
        sage: a.symmetries()
        no symmetry;  antisymmetry: (0, 1)
        sage: a[0,1] = 2
        sage: a[1,0]
        -2
        sage: a.comp()
        fully antisymmetric 2-indices components w.r.t. the coordinate frame (M, (d/dt,d/dx,d/dy,d/dz))
        sage: type(a.comp())
        <class 'sage.geometry.manifolds.component.CompFullyAntiSym'>

    Setting a component with repeated indices to a non-zero value results in an
    error::
    
        sage: a[1,1] = 3
        Traceback (most recent call last):
        ...
        ValueError: By antisymmetry, the component cannot have a nonzero value for the indices (1, 1)
        sage: a[1,1] = 0  # OK, albeit useless
        sage: a[1,2] = 3  # OK

    The expansion of a differential form with respect to a given coframe is 
    displayed via the method :meth:`view`::
    
        sage: a.view() # expansion with respect to the default coframe (dt, dx, dy, dz)
        a = 2 dt/\dx + 3 dx/\dy
        sage: latex(a.view()) # output for the notebook
        a = 2 \mathrm{d} t\wedge\mathrm{d} x + 3 \mathrm{d} x\wedge\mathrm{d} y

    Differential forms can be added or subtracted::
    
        sage: b = DiffForm(m, 2)
        sage: b[0,1], b[0,2], b[0,3] = (1,2,3)
        sage: s = a + b ; s
        2-form on the 4-dimensional manifold 'M'
        sage: a[:], b[:], s[:]
        (
        [ 0  2  0  0]  [ 0  1  2  3]  [ 0  3  2  3]
        [-2  0  3  0]  [-1  0  0  0]  [-3  0  3  0]
        [ 0 -3  0  0]  [-2  0  0  0]  [-2 -3  0  0]
        [ 0  0  0  0], [-3  0  0  0], [-3  0  0  0]
        )
        sage: s = a - b ; s
        2-form on the 4-dimensional manifold 'M'
        sage: s[:]
        [ 0  1 -2 -3]
        [-1  0  3  0]
        [ 2 -3  0  0]
        [ 3  0  0  0]

    An example of 3-form is the volume element on `\RR^3` in Cartesian 
    coordinates::
     
        sage: m = Manifold(3, 'R3', '\RR^3', start_index=1)                                
        sage: c_cart.<x,y,z> = m.chart('x y z')                                           
        sage: eps = DiffForm(m, 3, 'epsilon', r'\epsilon')
        sage: eps[1,2,3] = 1  # the only independent component
        sage: eps[:] # all the components are set from the previous line:
        [[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]]
        sage: eps.view()
        epsilon = dx/\dy/\dz
        
    Spherical components of the volume element from the tensorial 
    change-of-frame formula::

        sage: c_spher.<r,th,ph> = m.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: spher_to_cart = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: cart_to_spher = spher_to_cart.set_inverse(sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y, x))
        Check of the inverse coordinate transformation:
          r == r
          th == arctan2(r*sin(th), r*cos(th))
          ph == arctan2(r*sin(ph)*sin(th), r*cos(ph)*sin(th))
          x == x
          y == y
          z == z
        sage: eps.comp(c_spher.frame) # computation of the components in the spherical frame
        fully antisymmetric 3-indices components w.r.t. the coordinate frame (R3, (d/dr,d/dth,d/dph))
        sage: eps.comp(c_spher.frame)[1,2,3, c_spher]
        r^2*sin(th)
        sage: eps.view(c_spher.frame)
        epsilon = sqrt(x^2 + y^2 + z^2)*sqrt(x^2 + y^2) dr/\dth/\dph
        sage: eps.view(c_spher.frame, c_spher)
        epsilon = r^2*sin(th) dr/\dth/\dph
       
    The exterior product of two differential forms is performed via the method :meth:`wedge`::
    
        sage: a = OneForm(m, 'A')
        sage: a[:] = (x*y*z, -z*x, y*z)
        sage: b = OneForm(m, 'B')
        sage: b[:] = (cos(z), sin(x), cos(y))
        sage: ab = a.wedge(b) ; ab
        2-form 'A/\B' on the 3-dimensional manifold 'R3'
        sage: ab[:]
        [                         0  x*y*z*sin(x) + x*z*cos(z)  x*y*z*cos(y) - y*z*cos(z)]
        [-x*y*z*sin(x) - x*z*cos(z)                          0   -(x*cos(y) + y*sin(x))*z]
        [-x*y*z*cos(y) + y*z*cos(z)    (x*cos(y) + y*sin(x))*z                          0]

    The tensor product of a 1-form and a 2-form is not a 3-form but a tensor
    field of type (0,3) with less symmetries::
    
        sage: c = a*ab ; c 
        tensor field 'A*(A/\B)' of type (0,3) on the 3-dimensional manifold 'R3'
        sage: c.symmetries()  #  the antisymmetry is only w.r.t. the last two arguments:
        no symmetry;  antisymmetry: (1, 2)
        sage: d = ab*a ; d
        tensor field '(A/\B)*A' of type (0,3) on the 3-dimensional manifold 'R3'
        sage: d.symmetries()  #  the antisymmetry is only w.r.t. the first two arguments:
        no symmetry;  antisymmetry: (0, 1)

    The exterior derivative of a differential form is obtained by means of the 
    method :meth:`exterior_der`::
    
        sage: da = a.exterior_der() ; da
        2-form 'dA' on the 3-dimensional manifold 'R3'
        sage: da.view()
        dA = -(x + 1)*z dx/\dy - x*y dx/\dz + (x + z) dy/\dz
        sage: db = b.exterior_der() ; db
        2-form 'dB' on the 3-dimensional manifold 'R3'
        sage: db.view()
        dB = cos(x) dx/\dy + sin(z) dx/\dz - sin(y) dy/\dz
        sage: dab = ab.exterior_der() ; dab
        3-form 'd(A/\B)' on the 3-dimensional manifold 'R3'

    As a 3-form over a 3-dimensional manifold, d(A/\\B) is necessarily 
    proportional to the volume 3-form::
    
        sage: dab == dab[1,2,3]/eps[1,2,3]*eps
        True
        
    We may also check that the classical anti-derivation formula is fulfilled::
    
        sage: dab == da.wedge(b) - a.wedge(db)
        True
        
    The Lie derivative of a 2-form is a 2-form::
    
        sage: v = VectorField(m, 'v')             
        sage: v[:] = (y*z, -x*z, x*y)             
        sage: ab.lie_der(v)
        2-form on the 3-dimensional manifold 'R3'

    Let us check Cartan formula, which expresses the Lie derivative in terms
    of exterior derivatives::
    
        sage: ab.lie_der(v) == v.contract(0, ab.exterior_der(), 0) + (v.contract(0,ab,0)).exterior_der() 
        True
    
    """
    def __init__(self, domain, p, name=None, latex_name=None):
        TensorField.__init__(self, domain, 0, p, name, latex_name, 
                             antisym=range(p))
        DiffForm._init_derived(self) # initialization of derived quantities

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = str(self.rank) + "-form"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description

    def view(self, frame=None, chart=None):
        r"""
        Displays the differential form in terms of its expansion onto a given 
        coframe.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        INPUT:        
        
        - ``frame`` -- (default: None) vector frame with respect to which the 
          differential form is expanded; if none is provided, the domain's 
          default frame is assumed
        - ``chart`` -- (default: None) chart with respect to which the 
          components of the differential form in the selected frame are 
          expressed; if none is provided, the domain's default chart is assumed
          
        EXAMPLES:
        
        Display of a 2-form on `\RR^3`::
        
            sage: m = Manifold(3, 'R3', '\RR^3', start_index=1)
            sage: c_cart.<x,y,z> = m.chart('x y z')
            sage: a = DiffForm(m, 2, 'A')
            sage: a[1,2], a[2,3] = x*z, x^2+y^2
            sage: a.view() # expansion on the manifold's default coframe (dx, dy, dz)
            A = x*z dx/\dy + (x^2 + y^2) dy/\dz
            sage: latex(a.view()) # output for the notebook
            A = x z \mathrm{d} x\wedge\mathrm{d} y + \left( x^{2} + y^{2} \right) \mathrm{d} y\wedge\mathrm{d} z

        Display in a coframe different from the default one::
        
            sage: c_spher.<r,th,ph> = m.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi') # new coordinates 
            sage: spher_to_cart = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))  # the standard spherical coordinates
            sage: cart_to_spher = spher_to_cart.set_inverse(sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y, x), check=False)
            sage: a.view(c_spher.frame) # expansion on the coframe (dr, dth, dph) with the coefficients expressed in terms of the default coordinates (x,y,z)
            A = -sqrt(x^2 + y^2 + z^2)*sqrt(x^2 + y^2)*y dr/\dth + (x^3 + x*y^2 + x*z^2)*sqrt(x^2 + y^2) dth/\dph
            sage: a.view(c_spher.frame, c_spher) # expansion on the coframe (dr, dth, dph) with the coefficients expressed in terms of the spherical coordinates
            A = -r^3*sin(ph)*sin(th)^2 dr/\dth + r^4*cos(ph)*sin(th)^2 dth/\dph
            
        """
        from sage.misc.latex import latex
        from utilities import is_atomic, FormattedExpansion
        if frame is None:
            frame = self.domain.def_frame
        if chart is None:
            chart = self.domain.def_chart
        coframe = frame.coframe
        comp = self.comp(frame)
        terms_txt = []
        terms_latex = []
        n_con = self.tensor_type[0]
        for ind in comp.non_redundant_index_generator():
            coef = comp[[ind]].expr(chart)
            if coef != 0:
                bases_txt = []
                bases_latex = []
                for k in range(self.rank):
                    bases_txt.append(coframe[ind[k]].name)
                    bases_latex.append(latex(coframe[ind[k]]))
                basis_term_txt = "/\\".join(bases_txt)    
                basis_term_latex = r"\wedge".join(bases_latex)    
                if coef == 1:
                    terms_txt.append(basis_term_txt)
                    terms_latex.append(basis_term_latex)
                elif coef == -1:
                    terms_txt.append("-" + basis_term_txt)
                    terms_latex.append("-" + basis_term_latex)
                else:
                    coef_txt = repr(coef)
                    coef_latex = latex(coef)
                    if is_atomic(coef_txt):
                        terms_txt.append(coef_txt + " " + basis_term_txt)
                    else:
                        terms_txt.append("(" + coef_txt + ") " + 
                                         basis_term_txt)
                    if is_atomic(coef_latex):
                        terms_latex.append(coef_latex + basis_term_latex)
                    else:
                        terms_latex.append(r"\left(" + coef_latex + r"\right)" + 
                                           basis_term_latex)

        if terms_txt == []:
            expansion_txt = "0"
        else:
            expansion_txt = terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == "-":
                    expansion_txt += " - " + term[1:]
                else:
                    expansion_txt += " + " + term
        if terms_latex == []:
            expansion_latex = "0"
        else:
            expansion_latex = terms_latex[0]
            for term in terms_latex[1:]:
                if term[0] == "-":
                    expansion_latex += term
                else:
                    expansion_latex += "+" + term
        result = FormattedExpansion(self)            
        if self.name is None:
            result.txt = expansion_txt
        else:
            result.txt = self.name + " = " + expansion_txt
        if self.latex_name is None:
            result.latex = expansion_latex
        else:
            result.latex = latex(self) + " = " + expansion_latex
        return result


    def _new_comp(self, frame): 
        r"""
        Create some components in the given frame. 
                
        """
        from component import Components, CompFullyAntiSym
        if self.rank == 1: 
            return Components(frame, 1)
        else:
            return CompFullyAntiSym(frame, self.rank)

    def _new_instance(self):
        r"""
        Create a :class:`DiffForm` instance of the same degree. 
        
        """
        return DiffForm(self.domain, self.rank)

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        TensorField._init_derived(self)  
        self._exterior_derivative = None

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        TensorField._del_derived(self)
        self._exterior_derivative = None

    def degree(self):
        r"""
        Return the degree of the differential form. 
        """
        return self.rank

    def exterior_der(self, chart=None):
        r"""
        Compute the exterior derivative of the differential form. 
        
        INPUT:
        
        - ``chart`` -- (default: None): chart used for the computation; if none
          is provided, the computation is performed in a coordinate frame where
          the components are known, or computable by a component 
          transformation, privileging the domain's default chart.
        
        OUTPUT:
        
        - the exterior derivative of ``self``. 
        
        EXAMPLE:
        
        Exterior derivative of a 1-form on a 4-dimensional manifold::
        
            sage: m = Manifold(4, 'M')
            sage: c_txyz.<t,x,y,z> = m.chart('t x y z')           
            sage: a = OneForm(m, 'A')
            sage: a[:] = (t*x*y*z, z*y**2, x*z**2, x**2 + y**2)
            sage: da = a.exterior_der() ; da
            2-form 'dA' on the 4-dimensional manifold 'M'
            sage: da.view()
            dA = -t*y*z dt/\dx - t*x*z dt/\dy - t*x*y dt/\dz + (-2*y*z + z^2) dx/\dy + (-y^2 + 2*x) dx/\dz + (-2*x*z + 2*y) dy/\dz
            sage: latex(da)
            \mathrm{d}A
            
        The exterior derivative is nilpotent::
        
            sage: dda = da.exterior_der() ; dda
            3-form 'ddA' on the 4-dimensional manifold 'M'
            sage: dda.view()
            ddA = 0
            sage: dda == 0
            True

        """
        from sage.calculus.functional import diff
        from utilities import format_unop_txt, format_unop_latex
        if self._exterior_derivative is None:
            # A new computation is necessary:
            if chart is None:
                frame = self.pick_a_coord_frame()
                if frame is None:
                    raise ValueError("No coordinate frame could be found for " +
                                     "the differential form components.")
                chart = frame.chart
            else:
                frame = chart.frame
                if frame not in self.components:
                    raise ValueError("Components in the frame " + frame + 
                                     " have not been defined.")
            n = self.manifold.dim
            si = self.manifold.sindex
            sc = self.components[frame]
            dc = CompFullyAntiSym(frame, self.rank+1)
            for ind, val in sc._comp.items():
                for i in range(n):
                    ind_d = (i+si,) + ind
                    if len(ind_d) == len(set(ind_d)): # all indices are different
                          dc[[ind_d]] += val.function_chart(chart).diff(chart.xx[i])
            dc._del_zeros()
            # Name and LaTeX name of the result (rname and rlname):
            rname = format_unop_txt('d', self.name)
            rlname = format_unop_latex(r'\mathrm{d}', self.latex_name)
            # Final result
            self._exterior_derivative = DiffForm(self.domain, self.rank+1, 
                                                 rname, rlname)
            self._exterior_derivative.components[frame] = dc
        return self._exterior_derivative
 
 
    def wedge(self, other):
        r"""
        Compute the exterior product with another differential form. 
        
        INPUT:
        
        - ``other``: another differential form
        
        OUTPUT:
        
        - the exterior product self/\\other. 
        
        EXAMPLES:
        
        Exterior product of two 1-forms on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = m.chart('x y z')
            sage: a = OneForm(m, 'A')
            sage: a[:] = (x, y, z)
            sage: b = OneForm(m, 'B')
            sage: b[2] = z^2
            sage: a.view() ; b.view()
            A = x dx + y dy + z dz
            B = z^2 dz
            sage: h = a.wedge(b) ; h
            2-form 'A/\B' on the 3-dimensional manifold 'M'
            sage: h.view()
            A/\B = x*z^2 dx/\dz + y*z^2 dy/\dz
            sage: latex(h)
            A\wedge B
            sage: latex(h.view())
            A\wedge B = x z^{2} \mathrm{d} x\wedge\mathrm{d} z + y z^{2} \mathrm{d} y\wedge\mathrm{d} z

        The exterior product of two 1-forms is antisymmetric::
        
            sage: (b.wedge(a)).view()
            B/\A = -x*z^2 dx/\dz - y*z^2 dy/\dz
            sage: a.wedge(b) == - b.wedge(a)  
            True
            sage: (a.wedge(b) + b.wedge(a))[:]  # for the skeptical mind
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: a.wedge(a) == 0
            True

        The exterior product of a 2-form by a 1-form::
        
            sage: c = OneForm(m, 'C')
            sage: c[1] = y^3 ; c.view()
            C = y^3 dy
            sage: g = h.wedge(c) ; g
            3-form 'A/\B/\C' on the 3-dimensional manifold 'M'
            sage: g.view()
            A/\B/\C = -x*y^3*z^2 dx/\dy/\dz
            sage: g[:]
            [[[0, 0, 0], [0, 0, -x*y^3*z^2], [0, x*y^3*z^2, 0]], [[0, 0, x*y^3*z^2], [0, 0, 0], [-x*y^3*z^2, 0, 0]], [[0, -x*y^3*z^2, 0], [x*y^3*z^2, 0, 0], [0, 0, 0]]]

        The exterior product of a 2-form by a 1-form is symmetric::
        
            sage: h.wedge(c) == c.wedge(h)
            True

        """
        from utilities import is_atomic
        if not isinstance(other, DiffForm):
            raise TypeError("The second argument for the exterior product " + 
                            "must be a differential form.")
        if other.rank == 0:
            return other*self
        if self.rank == 0:
            return self*other
        frame = self.common_frame(other)
        if frame is None:
            raise ValueError("No common frame for the exterior product.")
        rank_r = self.rank + other.rank
        cmp_s = self.components[frame]
        cmp_o = other.components[frame]
        cmp_r = CompFullyAntiSym(frame, rank_r)
        for ind_s, val_s in cmp_s._comp.items():
            for ind_o, val_o in cmp_o._comp.items():
                ind_r = ind_s + ind_o
                if len(ind_r) == len(set(ind_r)): # all indices are different
                    cmp_r[ind_r] += val_s * val_o
        result = DiffForm(self.domain, rank_r)
        result.components[frame] = cmp_r
        if self.name is not None and other.name is not None:
            sname = self.name
            oname = other.name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            result.name = sname + '/\\' + oname
        if self.latex_name is not None and other.latex_name is not None:
            slname = self.latex_name
            olname = other.latex_name
            if not is_atomic(slname):
                slname = '(' + slname + ')'
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            result.latex_name = slname + r'\wedge ' + olname
        return result

    def hodge_star(self, metric):
        r"""
        Compute the Hodge dual of the differential form. 
        
        If ``self`` is a `p`-form `A`, its Hodge dual is the `(n-p)`-form
        `*A` defined by (`n` being the manifold's dimension)
        
        .. MATH::
            
            *A_{i_1\ldots i_{n-p}} = \frac{1}{p!} A_{k_1\ldots k_p}
                \epsilon^{k_1\ldots k_p}_{\qquad\ i_1\ldots i_{n-p}}
                
        where $\epsilon$ is the volume form associated with some 
        pseudo-Riemannian metric `g` on the manifold, and the indices 
        `k_1,\ldots, k_p` are raised with `g`. 
        
        INPUT:
        
        - ``metric``: the pseudo-Riemannian metric `g` defining the Hodge dual, 
          via the volume form `\epsilon`; must be an instance of :class:`Metric`
        
        OUTPUT:
        
        - the `(n-p)`-form `*A` 
        
        EXAMPLES:
        
        Hodge star of a 1-form in the Euclidean space `R^3`::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = m.chart('x y z')
            sage: g = Metric(m, 'g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: a = OneForm(m, 'A')
            sage: var('Ax Ay Az')
            (Ax, Ay, Az)
            sage: a[:] = (Ax, Ay, Az)
            sage: sa = a.hodge_star(g) ; sa
            2-form '*A' on the 3-dimensional manifold 'M'
            sage: sa.view()
            *A = Az dx/\dy - Ay dx/\dz + Ax dy/\dz
            sage: ssa = sa.hodge_star(g) ; ssa
            1-form '**A' on the 3-dimensional manifold 'M'
            sage: ssa.view()
            **A = Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Riemannian metric in dimension 3
            True
        
        Hodge star of a 0-form (scalar field) in `R^3`::
        
            sage: f = ScalarField(m, function('F',x,y,z), name='f')
            sage: sf = f.hodge_star(g) ; sf
            3-form '*f' on the 3-dimensional manifold 'M'
            sage: sf.view()
            *f = F(x, y, z) dx/\dy/\dz
            sage: ssf = sf.hodge_star(g) ; ssf
            scalar field '**f' on the 3-dimensional manifold 'M'
            sage: ssf.view()
            **f: (x, y, z) |--> F(x, y, z)
            sage: ssf == f # must hold for a Riemannian metric
            True
            
        Hodge star of a 0-form in Minkowksi spacetime::
        
            sage: m = Manifold(4, 'M')
            sage: X = m.chart('t x y z')
            sage: g = Metric(m, 'g', signature=2)
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: g.view()  # Minkowski metric
            g = -dt*dt + dx*dx + dy*dy + dz*dz
            sage: var('f0')
            f0
            sage: f = ScalarField(m, f0, name='f')
            sage: sf = f.hodge_star(g) ; sf 
            4-form '*f' on the 4-dimensional manifold 'M'
            sage: sf.view()
            *f = f0 dt/\dx/\dy/\dz
            sage: ssf = sf.hodge_star(g) ; ssf
            scalar field '**f' on the 4-dimensional manifold 'M'
            sage: ssf.view()
            **f: (t, x, y, z) |--> -f0
            sage: ssf == -f  # must hold for a Lorentzian metric             
            True

        Hodge star of a 1-form in Minkowksi spacetime::
        
            sage: a = OneForm(m, 'A')
            sage: var('At Ax Ay Az')
            (At, Ax, Ay, Az)
            sage: a[:] = (At, Ax, Ay, Az)
            sage: a.view()
            A = At dt + Ax dx + Ay dy + Az dz
            sage: sa = a.hodge_star(g) ; sa
            3-form '*A' on the 4-dimensional manifold 'M'
            sage: sa.view()
            *A = -Az dt/\dx/\dy + Ay dt/\dx/\dz - Ax dt/\dy/\dz - At dx/\dy/\dz
            sage: ssa = sa.hodge_star(g) ; ssa
            1-form '**A' on the 4-dimensional manifold 'M'
            sage: ssa.view()
            **A = At dt + Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Lorentzian metric in dimension 4
            True

        Hodge star of a 2-form in Minkowksi spacetime::
        
            sage: F = DiffForm(m, 2, 'F')    
            sage: var('Ex Ey Ez Bx By Bz')
            (Ex, Ey, Ez, Bx, By, Bz)
            sage: F[0,1], F[0,2], F[0,3] = -Ex, -Ey, -Ez
            sage: F[1,2], F[1,3], F[2,3] = Bz, -By, Bx
            sage: F[:]
            [  0 -Ex -Ey -Ez]
            [ Ex   0  Bz -By]
            [ Ey -Bz   0  Bx]
            [ Ez  By -Bx   0]
            sage: sF = F.hodge_star(g) ; sF
            2-form '*F' on the 4-dimensional manifold 'M'
            sage: sF[:]
            [  0  Bx  By  Bz]
            [-Bx   0  Ez -Ey]
            [-By -Ez   0  Ex]
            [-Bz  Ey -Ex   0]
            sage: ssF = sF.hodge_star(g) ; ssF
            2-form '**F' on the 4-dimensional manifold 'M'
            sage: ssF[:]   
            [  0  Ex  Ey  Ez]
            [-Ex   0 -Bz  By]
            [-Ey  Bz   0 -Bx]
            [-Ez -By  Bx   0]
            sage: ssF.view()
            **F = Ex dt/\dx + Ey dt/\dy + Ez dt/\dz - Bz dx/\dy + By dx/\dz - Bx dy/\dz
            sage: F.view()
            F = -Ex dt/\dx - Ey dt/\dy - Ez dt/\dz + Bz dx/\dy - By dx/\dz + Bx dy/\dz
            sage: ssF == -F  # must hold for a Lorentzian metric in dimension 4
            True

        Test of the standard identity
        
        .. MATH::
            
            *(A\wedge B) = \epsilon(A^\sharp, B^\sharp, ., .)
            
        where `A` and `B` are any 1-forms and `A^\sharp` and `B^\sharp` the 
        vectors associated to them by the metric `g` (index raising)::

            sage: b = OneForm(m, 'B')
            sage: var('Bt Bx By Bz')
            (Bt, Bx, By, Bz)
            sage: b[:] = (Bt, Bx, By, Bz) ; b.view()
            B = Bt dt + Bx dx + By dy + Bz dz
            sage: epsilon = g.volume_form()
            sage: (a.wedge(b)).hodge_star(g) == epsilon.contract(0, a.up(g), 0).contract(0, b.up(g), 0)
            True

        """
        from sage.functions.other import factorial
        from utilities import format_unop_txt, format_unop_latex
        p = self.rank
        eps = metric.volume_form(p)
        if p == 0:
            resu = self * eps
        else:
            resu = self.contract(0, eps, 0)
            for j in range(1, p):
                resu = resu.self_contract(0, p-j)
            if p > 1:
                resu = resu / factorial(p)
        # Name and LaTeX name of the result:
        resu.name = format_unop_txt('*', self.name)
        resu.latex_name = format_unop_latex(r'\star ', self.latex_name)
        return resu
        
        
#******************************************************************************

class OneForm(DiffForm):
    r"""
    Class for 1-forms on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the 1-form is defined
    - ``name`` -- (default: None) name given to the 1-form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 1-form; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    A 1-form on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')                      
        sage: c_xyz.<x,y,z> = m.chart('x y z')
        sage: om = OneForm(m, 'omega', r'\omega') ; om  
        1-form 'omega' on the 3-dimensional manifold 'M'

    A 1-form is of course a differential form::
    
        sage: isinstance(om, DiffForm)
        True
        sage: om.tensor_type
        (0, 1)
        
    Setting the components w.r.t. the manifold's default frame::
    
        sage: om[:] = (2*z, x, x-y)
        sage: om[:]
        [2*z, x, x - y]
        sage: om.view()
        omega = 2*z dx + x dy + (x - y) dz
        
    A 1-form acts on vector fields::
    
        sage: v = VectorField(m, 'V')
        sage: v[:] = (x, 2*y, 3*z)
        sage: om(v)
        scalar field 'omega(V)' on the 3-dimensional manifold 'M'
        sage: om(v).view()
        omega(V): (x, y, z) |--> 2*x*y + (5*x - 3*y)*z
        sage: latex(om(v))
        \omega\left(V\right)

    The tensor product of two 1-forms is a tensor field of type (0,2)::
    
        sage: a = OneForm(m, 'A')                       
        sage: a[:] = (1, 2, 3)                          
        sage: b = OneForm(m, 'B')
        sage: b[:] = (6, 5, 4)
        sage: c = a*b ; c       
        tensor field 'A*B' of type (0,2) on the 3-dimensional manifold 'M'
        sage: c[:]
        [ 6  5  4]
        [12 10  8]
        [18 15 12]
        sage: c.symmetries()    # c has no symmetries:
        no symmetry;  no antisymmetry
        
    The exterior product of two 1-forms is a 2-form::
    
        sage: d = a.wedge(b) ; d
        2-form 'A/\B' on the 3-dimensional manifold 'M'
        sage: d[:]
        [  0  -7 -14]
        [  7   0  -7]
        [ 14   7   0]
        sage: d.symmetries()
        no symmetry;  antisymmetry: (0, 1)

    We can check the standard formula relating the exterior product to the
    tensor product::
    
        sage: a.wedge(b) == a*b - b*a
        True

    """
    def __init__(self, domain, name=None, latex_name=None):
        DiffForm.__init__(self, domain, 1, name, latex_name)
        
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "1-form"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description

    def _new_comp(self, frame): 
        r"""
        Create some components in the given frame. 
                
        """
        from component import Components
        return Components(frame, 1)

    def _new_instance(self):
        r"""
        Create a :class:`OneForm` instance.
        
        """
        return OneForm(self.domain)

    def __call__(self, vector):
        r"""
        The 1-form acting on a vector field.
        
        INPUT:
        
        - ``vector`` -- a vector field `v` (instance of :class:`VectorField`)
        
        OUTPUT:
        
        - the scalar field `\langle \omega, v \rangle`, as an instance of 
          :class:`ScalarField`
          
        """
        from vectorfield import VectorField
        if not isinstance(vector, VectorField):
            raise TypeError("The argument must be a vector field.")
        frame = self.common_frame(vector)
        if frame is None:
            raise ValueError("No common frame for the components.")
        omega = self.components[frame]
        vv = vector.components[frame]
        si = self.manifold.sindex
        resu = omega[[si]]*vv[[si]]
        for i in self.manifold.irange(si+1):
            resu += omega[[i]]*vv[[i]]
        # Name of the output:
        resu.name = None
        if self.name is not None and vector.name is not None:
            resu.name = self.name + "(" + vector.name + ")"
        # LaTeX symbol for the output:
        resu.latex_name = None
        if self.latex_name is not None and vector.latex_name is not None:
            resu.latex_name = self.latex_name + r"\left(" + \
                              vector.latex_name + r"\right)"
        return resu


        
        
        
        
        
        
        
        
        
        
        
        
        
        
