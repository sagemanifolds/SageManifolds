r"""
Tensor fields

The class :class:`TensorField` implements tensor fields on differentiable 
manifolds over `\RR`. 

A tensor field of type `(k,\ell)` is a field of multilinear maps:

.. MATH::

    \underbrace{T_p^*(M)\times\cdots\times T_p^*(M)}_{k\ \; \mbox{times}}
    \times \underbrace{T_p(M)\times\cdots\times T_p(M)}_{\ell\ \; \mbox{times}}
    \longrightarrow \RR
    
where `T_p(M)` stands for the tangent space at the point `p` on the
manifold `M` and `T_p^*(M)` for its dual vector space. The integer `k+\ell`
is called the tensor rank. 

Various derived classes of :class:`TensorField` are devoted to specific tensor
fields:

* :class:`ScalarField` for scalar fields (rank-0 tensor fields)
* :class:`VectorField` for vector fields (rank-1 contravariant tensor fields)
* :class:`OneForm` for 1-forms (rank-1 convariant tensor fields)
* :class:`EndomorphismField` for fields of endomorphisms (type (1,1) tensor 
  fields)
* :class:`SymBilinFormField` for fields of symmetric bilinear forms (rank-2
  symmetric covariant tensor fields)
* :class:`DiffForm` for differential forms (fully antisymmetric covariant 
  tensor fields)


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

EXAMPLES:

    A tensor field of type (1,1) on a 2-dimensional manifold::

        sage: m = Manifold(2, 'M', start_index=1)
        sage: c_xy = Chart(m, 'x y', 'coord_xy')
        sage: t = TensorField(m, 1, 1, 'T') ; t
        tensor field 'T' of type (1,1) on the 2-dimensional manifold 'M'
        sage: t.rank
        2

    A just-created tensor field has no components::

        sage: t.components
        {}

    Components w.r.t. the manifold's default frame are created by providing the
    relevant indices inside square brackets::

        sage: t[1,1] = x^2

    Unset components are initialized to zero::

        sage: t[:]  # list of components w.r.t. the manifold's default vector frame
        [x^2   0]
        [  0   0]

    The full set of components w.r.t. a given vector frame is returned by the 
    method :meth:`comp`; it is an instance of the class :class:`Components`::
    
        sage: t.comp('coord_xy_b')  
        2-indices components w.r.t. the coordinate basis 'coord_xy_b' (d/dx,d/dy)
        sage: print type(t.comp('coord_xy_b'))
        <class 'sage.geometry.manifolds.component.Components'>

    The name of the vector frame can be skipped, it is then assumed to be the
    manifold's default frame::
    
        sage: m.default_frame()
        coordinate basis 'coord_xy_b' (d/dx,d/dy)
        sage: t.comp() is t.comp('coord_xy_b')
        True

    Individual components w.r.t. the manifold's default frame are accessed by 
    listing their indices inside double square brackets; they are scalar
    fields on the manifold, and therefore instances of the class 
    :class:`ScalarField`::

        sage: t[[1,1]]
        scalar field on the 2-dimensional manifold 'M'
        sage: t[[1,1]].expr()
        x^2
        sage: t[[1,2]]
        zero scalar field on the 2-dimensional manifold 'M'
        sage: t[[1,2]].expr()
        0
        
    A direct access to the coordinate expression of some component is obtained
    via the single square brackets::
    
        sage: t[1,1] 
        x^2
        sage: t[1,1] is t[[1,1]].function_chart() # the coordinate function
        True
        sage: t[1,1] is t[[1,1]].function_chart('coord_xy')
        True
        sage: t[1,1] == t[[1,1]].expr() # check the value of the coordinate function
        True
        sage: t[1,1].expr() is t[[1,1]].expr() # the symbolic expression
        True

    In other words, the single square brackets return an instance of 
    :class:`FunctionChart` that is the coordinate function representing the 
    component in some chart (by default, the manifold's default chart)::
    
        sage: print type(t[1,1])    # single bracket --> FunctionChart
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: print type(t[[1,1]])  # double bracket --> ScalarField
        <class 'sage.geometry.manifolds.scalarfield.ScalarField'>
    
    Expressions in a chart different from the manifold's default one are 
    obtained by specifying the chart name as the last argument inside the
    single square brackets::
    
        sage: c_uv = Chart(m, 'u v', 'coord_uv')
        sage: xy_to_uv = CoordChange(c_xy, c_uv, x+y, x-y)  
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: t[1,1, 'coord_uv'] 
        1/4*u^2 + 1/2*u*v + 1/4*v^2

    Note that ``t[1,1, 'coord_uv']`` is the component of the tensor t w.r.t. to 
    the vector frame 'coord_xy_b' (the coordinate basis associated to the 
    chart (x,y)) expressed in terms of the coordinates (u,v). Indeed, 
    ``t[1,1, 'coord_uv']`` is a shortcut for 
    ``t.comp('coord_xy_b')[[1,1]].function_chart('coord_uv')``::
        
        sage: t[1,1, 'coord_uv'] is t.comp('coord_xy_b')[[1,1]].function_chart('coord_uv')
        True
    
    Similarly, ``t[1,1]`` is a shortcut for 
    ``t.comp('coord_xy_b')[[1,1]].function_chart('coord_xy')``::
    
        sage: t[1,1] is t.comp('coord_xy_b')[[1,1]].function_chart('coord_xy')            
        True
        sage: t[1,1] is t.comp()[[1,1]].function_chart()  # since 'coord_xy_b' and 'coord_xy' are the manifold's default values                
        True

    Internally, the components are stored as a dictionary (attribute 
    :attr:`_comp` of the class :class:`Components`) whose
    keys are the indices. Only the non-zero components and non-redundant
    components (in case of symmetries) are stored::

        sage: t.comp()._comp
        {(1, 1): scalar field on the 2-dimensional manifold 'M'}

    All the components can be set at once via [:]::
    
        sage: t[:] = [[1, -x], [x*y, 2]]
        sage: t[:]
        [  1  -x]
        [x*y   2]
        
    The different sets of components, corresponding to representations of the
    tensor in different vector frames, are stored in the dictionary 
    :attr:`components`, each item being an instance of the class 
    :class:`Components`::
    
        sage: t.components
        {'coord_xy_b': 2-indices components w.r.t. the coordinate basis 'coord_xy_b' (d/dx,d/dy)}
        sage: print type(t.components['coord_xy_b'])
        <class 'sage.geometry.manifolds.component.Components'>
        sage: print type(t.comp())
        <class 'sage.geometry.manifolds.component.Components'>
        sage: t.comp() is t.components['coord_xy_b']
        True

    To set the components in a vector frame different from the manifold's 
    default one, the method :meth:`set_comp` must be employed::

        sage: e = VectorFrame(m, 'e')
        sage: t.set_comp('e')[1,1], t.set_comp('e')[1,2] = (x+y, 0)
        sage: t.set_comp('e')[2,1], t.set_comp('e')[2,2] = (y, -3*x)
        sage: t.comp('e')
        2-indices components w.r.t. the vector frame 'e'
        sage: t.comp('e')[:]
        [x + y     0]
        [    y  -3*x]

    All the components in some frame can be set at once, via the operator
    [:]::

        sage: t.set_comp('e')[:] = [[x+y, 0], [y, -3*x]]
        sage: t.comp('e')[:]  # same as above:
        [x + y     0]
        [    y  -3*x]
    
    To avoid any insconstency between the various components, the default
    behaviour of :meth:`set_comp` is to delete the components in other frames. 
    Accordingly, the components in the frame 'coord_xy_b' have been deleted::
    
        sage: t.components
        {'e': 2-indices components w.r.t. the vector frame 'e'}

    To keep the other components, the keyword ``delete_others`` has to be set to 
    False::
    
        sage: t = TensorField(m, 1, 1, 'T')  # Let us restart 
        sage: t[:] = [[1, -x], [x*y, 2]]  # by first setting the components in the frame 'coord_xy_b'
        sage: # We now set the components in the frame 'e' with delete_others=False:
        sage: t.set_comp('e', delete_others=False)[:] = [[x+y, 0], [y, -3*x]]
        sage: t.components  # Both set of components are present:
        {'e': 2-indices components w.r.t. the vector frame 'e', 'coord_xy_b': 2-indices components w.r.t. the coordinate basis 'coord_xy_b' (d/dx,d/dy)}

    The expansion of the tensor field in a given frame is displayed via the 
    method :meth:`show` (the symbol * stands for tensor product)::
    
        sage: t.show()  # expansion in the manifold's default frame
        T = d/dx*dx - x d/dx*dy + x*y d/dy*dx + 2 d/dy*dy
        sage: t.show('e')
        T = (x + y) e_1*e^1 + y e_2*e^1 - 3*x e_2*e^2

    A tensor field acts as a multilinear map on 1-forms and vector fields; 
    in the present case, T being of type (1,1), it acts on pairs 
    (1-form, vector)::
    
        sage: a = OneForm(m, 'a')
        sage: a[:] = (1, x)
        sage: v = VectorField(m, 'V')
        sage: v[:] = (y, 2)
        sage: t(a,v)
        scalar field 'T(a,V)' on the 2-dimensional manifold 'M'
        sage: t(a,v).expr()
        x^2*y^2 + 2*x + y
        sage: latex(t(a,v))
        T\left(a,V\right)
    
    Check by means of the component expression of t(a,v)::
    
        sage: t[1,1]*a[1]*v[1] + t[1,2]*a[1]*v[2] + t[2,1]*a[2]*v[1] + t[2,2]*a[2]*v[2] - t(a,v).expr()
        0

    A scalar field (rank-0 tensor field)::
    
        sage: f = ScalarField(m, x*y + 2, name='f') ; f 
        scalar field 'f' on the 2-dimensional manifold 'M'
        sage: isinstance(f, TensorField)
        True
        sage: f.tensor_type
        (0, 0)
        
    As differential mappings from the manifold to `\RR`, scalar fields also 
    inherit from the class :class:`DiffMapping`::
    
        sage: isinstance(f, DiffMapping)
        True
        sage: f.manifold1 # the start manifold
        2-dimensional manifold 'M'
        sage: f.manifold2 # the target manifold
        field R of real numbers

    They act on points on the manifold (as any instance of 
    :class:`DiffMapping`)::
    
        sage: p = Point(m, (1,2))
        sage: f(p)
        4
        
    A vector field (rank-1 contravariant tensor field)::
    
        sage: v = VectorField(m, 'v') ; v
        vector field 'v' on the 2-dimensional manifold 'M'
        sage: v.tensor_type
        (1, 0)
        sage: v[1], v[2] = -x, y
        sage: v.show()
        v = -x d/dx + y d/dy        

    A field of symmetric bilinear forms::
    
        sage: q = SymBilinFormField(m, 'Q') ; q
        field of symmetric bilinear forms 'Q' on the 2-dimensional manifold 'M'
        sage: q.tensor_type
        (0, 2)

    The components of a symmetric bilinear form are dealt by the subclass 
    :class:`CompFullySym` of the class :class:`Components`, which takes into 
    account the symmetry between the two indices::
    
        sage: q[1,1], q[1,2], q[2,2] = (0, -x, y) # no need to set the component (2,1)
        sage: print type(q.comp())
        <class 'sage.geometry.manifolds.component.CompFullySym'>
        sage: q[:] # note that the component (2,1) is equal to the component (1,2)
        [ 0 -x]
        [-x  y]
        sage: q.show()
        Q = -x dx*dy - x dy*dx + y dy*dy
    
    Internally (dictionary :attr:`_comp` of the class :class:`Components`), only
    the non-zero and non-redundant components are stored::
    
        sage: q.comp()._comp
        {(1, 2): scalar field on the 2-dimensional manifold 'M',
        (2, 2): scalar field on the 2-dimensional manifold 'M'}
        sage: q.comp()._comp[(1,2)].expr()
        -x
        sage: q.comp()._comp[(2,2)].expr()
        y

    More generally, tensor symmetries or antisymmetries can be specified via
    the keywords ``sym`` and ``antisym``. For instance a rank-4 covariant 
    tensor symmetric with respect to its first two arguments and 
    antisymmetric with respect to its last two ones is declared as follows::
    
        sage: t = TensorField(m, 0, 4, 'T', sym=(0,1), antisym=(2,3))
        sage: t[1,2,1,2] = 3
        sage: t[2,1,1,2] # check of the symmetry with respect to the first 2 indices
        3
        sage: t[1,2,2,1] # check of the antisymmetry with respect to the last 2 indices
        -3

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
from sage.rings.integer import Integer
from manifold import Manifold
from component import Components, CompWithSym, CompFullySym, CompFullyAntiSym

class TensorField(SageObject):
    r"""
    Base class for tensor fields on a differentiable manifold.
    
    A tensor field of type `(k,\ell)` is a field of multilinear maps:

    .. MATH::

        \underbrace{T_p^*(M)\times\cdots\times T_p^*(M)}_{k\ \; \mbox{times}}
        \times \underbrace{T_p(M)\times\cdots\times T_p(M)}_{\ell\ \; \mbox{times}}
        \longrightarrow \RR
    
    where `T_p(M)` stands for the tangent space at the point `p` on the
    manifold `M` and `T_p^*(M)` for its dual vector space. The integer `k+\ell`
    is called the tensor rank. 
    
    INPUT:
    
    - ``manifold`` -- the manifold on which the tensor field is defined
    - ``k`` -- the contravariant rank, the tensor type being (k,l)
    - ``l`` -- the covariant rank, the tensor type being (k,l)
    - ``name`` -- (default: None) name given to the tensor field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
      if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      tensor arguments: each symmetry is described by a tuple containing 
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:
        * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
        * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
          arguments and a symmetry between the 2nd, 4th and 5th arguments.
    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the arguments, with the same convention as for ``sym``. 


    EXAMPLES:

    A tensor field of type (2,0) on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = Chart(m, 'x y z', 'xyz-coord')
        sage: t = TensorField(m, 2, 0, 'T') ; t
        tensor field 'T' of type (2,0) on the 3-dimensional manifold 'M'

    The components with respect to the manifold's default frame are set or read
    by means of square brackets::
    
        sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
        sage: for i in range(3):
        ...       for j in range(3):
        ...           t[i,j] = (i+1)**(j+1)
        ...
        sage: [[ t[i,j] for j in range(3)] for i in range(3)]
        [[1, 1, 1], [2, 4, 8], [3, 9, 27]]
    
    A shortcut for the above is using [:]::
    
        sage: t[:]
        [ 1  1  1]
        [ 2  4  8]
        [ 3  9 27]

    The components with respect to another frame are set via the method
    :meth:`set_comp` and read via the method :meth:`comp`; both return an 
    instance of :class:`Components`::
    
        sage: f = VectorFrame(m, 'f')  # a new frame defined on M, in addition to e
        sage: t.set_comp('f')[0,0] = -3
        sage: t.comp('f')
        2-indices components w.r.t. the vector frame 'f'
        sage: t.comp('f')[0,0]
        -3
        sage: t.comp('f')[:]  # the full list of components
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]

    To avoid any insconstency between the various components, the default
    behaviour of :meth:`set_comp` is to delete the components in other frames. 
    Accordingly, the components in the frame e have been deleted::
    
        sage: t.components
        {'f': 2-indices components w.r.t. the vector frame 'f'}

    To keep the other components, one must set the keyword ``delete_others`` to 
    False::
    
        sage: t = TensorField(m, 2, 0, 'T')  # Let us restart
        sage: t[0,0] = 2                     # sets the components in the frame e
        sage: # We now set the components in the frame f with delete_others=False:
        sage: t.set_comp('f', delete_others=False)[0,0] = -3
        sage: # The components w.r.t. the frame e have been kept: 
        sage: t.components
        {'e': 2-indices components w.r.t. the vector frame 'e', 'f': 2-indices components w.r.t. the vector frame 'f'}

    The basic attributes of :class:`TensorField` are :attr:`manifold`, rank 
    (rank), :attr:`tensor_type` (the pair (k,l)) and :attr:`components` (the
    dictionary of the components w.r.t. various frames)::

        sage: t.manifold
        3-dimensional manifold 'M'
        sage: t.rank
        2
        sage: t.tensor_type
        (2, 0)
        sage: t.components
        {'e': 2-indices components w.r.t. the vector frame 'e', 'f': 2-indices components w.r.t. the vector frame 'f'}

    Symmetries and antisymmetries are declared via the keywords ``sym`` and
    ``antisym``. For instance, a rank-6 covariant tensor that is symmetric with
    respect to its 1st and 3rd arguments and antisymmetric with respect to the 
    2nd, 5th and 6th arguments is set up as follows::
    
        sage: a = TensorField(m, 0, 6, 'T', sym=(0,2), antisym=(1,4,5))
        sage: a[0,0,1,0,1,2] = 3
        sage: a[1,0,0,0,1,2] # check of the symmetry
        3
        sage: a[0,1,1,0,0,2], a[0,1,1,0,2,0] # check of the antisymmetry
        (-3, 3)  

    Multiple symmetries or antisymmetries are allowed; they must then be 
    declared as a list. For instance, a rank-4 covariant tensor that is 
    antisymmetric with respect to its 1st and 2nd arguments and with respect to
    its 3rd and 4th argument must be declared as::
    
        sage: r = TensorField(m, 0, 4, 'T', antisym=[(0,1), (2,3)])
        sage: r[0,1,2,0] = 3
        sage: r[1,0,2,0] # first antisymmetry
        -3
        sage: r[0,1,0,2] # second antisymmetry
        -3
        sage: r[1,0,0,2] # both antisymmetries acting
        3
    
    Tensor fields of the same type can be added and subtracted::
    
        sage: a = TensorField(m, 2, 0)
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: b = TensorField(m, 2, 0)
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + 2*b ; s
        tensor field of type (2,0) on the 3-dimensional manifold 'M'
        sage: a[:], (2*b)[:], s[:]
        (
        [1 2 3]  [ 8  0 14]  [ 9  2 17]
        [0 0 0]  [ 0 10  0]  [ 0 10  0]
        [0 0 0], [ 0  0 12], [ 0  0 12]
        )
        sage: s = a - b ; s
        tensor field of type (2,0) on the 3-dimensional manifold 'M'
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [-3  2 -4]
        [0 0 0]  [0 5 0]  [ 0 -5  0]
        [0 0 0], [0 0 6], [ 0  0 -6]
        )

    Symmetries are preserved by the addition whenever it is possible::
    
        sage: a = TensorField(m, 2, 0, sym=(0,1))
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [0 0 6], [ 3  0  6]
        )
        sage: a.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: b.symmetries()
        no symmetry;  no antisymmetry
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: # let us now make b symmetric:
        sage: b = TensorField(m, 2, 0, sym=(0,1))
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [7 0 6], [10  0  6]
        )
        sage: s.symmetries()  # s is symmetric because both a and b are
        symmetry: (0, 1);  no antisymmetry

    The tensor product is taken with the operator \*::
    
        sage: c = a*b ; c
        tensor field of type (4,0) on the 3-dimensional manifold 'M'
        sage: c.symmetries()  # since a and b are both symmetric, a*b has two symmetries:
        symmetries: [(0, 1), (2, 3)];  no antisymmetry

    The tensor product of two fully contravariant tensors is not symmetric in 
    general::
    
        sage: a*b == b*a
        False

    The tensor product of a fully contravariant tensor by a fully covariant one
    is symmetric::
    
        sage: d = DiffForm(m, 2)  # a fully covariant tensor field
        sage: d[0,1], d[0,2], d[1,2] = (3, 2, 1)
        sage: s = a*d ; s 
        tensor field of type (2,2) on the 3-dimensional manifold 'M'
        sage: s.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: s1 = d*a ; s1 
        tensor field of type (2,2) on the 3-dimensional manifold 'M'
        sage: s1.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: d*a == a*d
        True

    """
    def __init__(self, manifold, k, l, name=None, latex_name=None, sym=None, 
                 antisym=None):
        if not isinstance(manifold, Manifold):
            raise TypeError("The first argument must be a manifold.")
        self.manifold = manifold
        self.tensor_type = (k,l)
        self.rank = k+l
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.components = {}    # components not set yet
        # Treatment of symmetries declarations:
        self.sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):  
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                         " not in [0," + str(self.rank-1) + "]")
                    self.sym.append(tuple(isym))       
        self.antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):  
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                         " not in [0," + str(self.rank-1) + "]")
                    self.antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self.sym:
            index_list += isym
        for isym in self.antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("Incompatible lists of symmetries: the same " + 
                             "position appears more then once.")
        # Initialization of derived quantities:
        self._init_derived() 

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "tensor field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " of type (%s,%s)" % (str(self.tensor_type[0]), 
                                             str(self.tensor_type[1]))
        description += " on the " + str(self.manifold)
        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{no symbol}'
        else:
           return self.latex_name

    def set_name(self, name, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of the tensor.

        INPUT:
        
        - ``name`` -- name given to the tensor field
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor 
        field; if none is provided, the LaTeX symbol is set to ``name``

        """
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
       

    def show(self, frame_name=None, chart_name=None):
        r"""
        Displays the tensor field in terms of its expansion onto a given frame.
        
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        INPUT:
                
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame with respect to which the tensor field is expanded; if 
          none is provided, the manifold's default frame is assumed
        - ``chart_name`` -- (default: None) string containing the name of the 
          chart with respect to which the components of the tensor field in the 
          selected frame are expressed; if none is provided, the manifold's 
          default chart is assumed

        EXAMPLES:
        
        Displaying expansions on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: v = VectorField(m, 'v')
            sage: v[1], v[2] = y+1, -x
            sage: v.show()  # expansion in the manifold's default frame
            v = (y + 1) d/dx - x d/dy
            sage: latex(v.show()) # display in the notebook
            v = \left( y + 1 \right) \frac{\partial}{\partial x } -x \frac{\partial}{\partial y }
            sage: e = VectorFrame(m, 'e')
            sage: v.set_comp('e', delete_others=False)[:] = 2, x*y
            sage: v.show('e')
            v = 2 e_1 + x*y e_2
            sage: latex(v.show('e'))
            v = 2 e_1 + x y e_2
            sage: t = TensorField(m, 1, 1, 'T')
            sage: t[:] = [[1, -x], [x*y, 2]]
            sage: t.set_comp('e', delete_others=False)[:] = [[x+y, 0], [y, -3*x]]
            sage: t.show() # the symbol * stands for the tensor product
            T = d/dx*dx - x d/dx*dy + x*y d/dy*dx + 2 d/dy*dy
            sage: t.show('e')
            T = (x + y) e_1*e^1 + y e_2*e^1 - 3*x e_2*e^2
            sage: latex(t.show('e'))
            T = \left( x + y \right) e_1\otimes e^1 + y e_2\otimes e^1 -3 \, x e_2\otimes e^2
            sage: q = SymBilinFormField(m, 'Q')
            sage: q[1,1], q[1,2], q[2,2] = (0, -x, y)
            sage: q.show()
            Q = -x dx*dy - x dy*dx + y dy*dy
            sage: latex(q.show())
            Q = -x \mathrm{d} x\otimes \mathrm{d} y -x \mathrm{d} y\otimes \mathrm{d} x + y \mathrm{d} y\otimes \mathrm{d} y
            
        """
        from sage.misc.latex import latex
        from utilities import is_atomic, FormattedExpansion
        if frame_name is None:
            frame_name = self.manifold.def_frame.name
        if chart_name is None:
            chart_name = self.manifold.def_chart.name
        frame = self.manifold.frames[frame_name]
        coframe = frame.coframe
        comp = self.comp(frame_name)
        terms_txt = []
        terms_latex = []
        n_con = self.tensor_type[0]
        for ind in self.manifold.index_generator(self.rank):
            coef = comp[[ind]].expr(chart_name)
            if coef != 0:
                bases_txt = []
                bases_latex = []
                for k in range(n_con):
                    bases_txt.append(frame(ind[k]).name)
                    bases_latex.append(latex(frame(ind[k])))
                for k in range(n_con, self.rank):
                    bases_txt.append(coframe(ind[k]).name)
                    bases_latex.append(latex(coframe(ind[k])))
                basis_term_txt = "*".join(bases_txt)    
                basis_term_latex = r"\otimes ".join(bases_latex)    
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
    
    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries.
        
        EXAMPLES:
        
        Various symmetries / antisymmetries for a rank-4 tensor::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = Chart(m, 'x y z', 'xyz-coord')
            sage: t = TensorField(m, 4, 0, 'T') # no symmetry declared
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: t = TensorField(m, 4, 0, 'T', sym=(0,1))
            sage: t.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: t = TensorField(m, 4, 0, 'T', sym=[(0,1), (2,3)])
            sage: t.symmetries()
            symmetries: [(0, 1), (2, 3)];  no antisymmetry
            sage: t = TensorField(m, 4, 0, 'T', sym=(0,1), antisym=(2,3))
            sage: t.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            
        """
        if len(self.sym) == 0:
            s = "no symmetry; "
        elif len(self.sym) == 1:
            s = "symmetry: " + str(self.sym[0]) + "; "
        else:
            s = "symmetries: " + str(self.sym) + "; " 
        if len(self.antisym) == 0:
            a = "no antisymmetry"
        elif len(self.antisym) == 1:
            a = "antisymmetry: " + str(self.antisym[0])
        else:
            a = "antisymmetries: " + str(self.antisym)   
        print s, a
         
    def _new_comp(self, frame_name): 
        r"""
        Create some components in the given frame. 
        
        This method, to be called by :meth:`comp`, must be redefined by derived 
        classes to adapt the output to the relevant subclass of 
        :class:`Components`.
        
        """
        if self.sym == [] and self.antisym == []:
            return Components(self.manifold, self.rank, frame_name)
        for isym in self.sym:
            if len(isym) == self.rank:
                return CompFullySym(self.manifold, self.rank, frame_name)
        for isym in self.antisym:
            if len(isym) == self.rank:
                return CompFullyAntiSym(self.manifold, self.rank, frame_name)
        return CompWithSym(self.manifold, self.rank, frame_name, sym=self.sym,
                           antisym=self.antisym)        
    
    def _new_instance(self):
        r"""
        Create a :class:`TensorField` instance of the same tensor type and 
        with the same symmetries.

        This method must be redefined by derived classes of 
        :class:`TensorField`.
        
        """
        return TensorField(self.manifold, self.tensor_type[0], 
                           self.tensor_type[1], sym=self.sym, 
                           antisym=self.antisym)
        
    def _init_derived(self):
        r"""
        Initialises the derived quantities
        """
        self._lie_derivatives = {} # collection of Lie derivatives of self

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        # First deletes any reference to self in the vectors' dictionary:
        for vid, val in self._lie_derivatives.items():
            del val[0]._lie_der_along_self[id(self)]
        self._lie_derivatives.clear()

    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The name and the derived quantities are not copied. 
        
        EXAMPLES:

        Copy on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')  
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: t = TensorField(m, 1, 1, 'T') ; t
            tensor field 'T' of type (1,1) on the 2-dimensional manifold 'M'
            sage: t[:] = [[1, x+y], [x*y, 2]]
            sage: t1 = t.copy() ; t1  # note that the name 'T' is not copied
            tensor field of type (1,1) on the 2-dimensional manifold 'M'
            sage: t1[:]
            [    1 x + y]
            [  x*y     2] 
            sage: t1 == t
            True
            sage: t1 is t  # the copy is not the same object
            False
            sage: t1.comp() is t.comp()  # the components are also different objects
            False

        The precise tensor kind is preserved, and thus the symmetries::
        
            sage: t = SymBilinFormField(m)
            sage: t[0,0], t[0,1], t[1,1] = (x, 2, x*y)
            sage: t1 = t.copy() ; t1
            field of symmetric bilinear forms on the 2-dimensional manifold 'M'
            sage: t1[:]
            [  x   2]
            [  2 x*y]
            sage: t1 == t
            True

        """
        result = self._new_instance()
        for frame_name, comp in self.components.items():
             result.components[frame_name] = comp.copy()
        return result


    def comp(self, frame_name=None, from_frame=None):
        r"""
        Return the components in a given frame.
        
        If the components are not known already, they are computed by the tensor
        change-of-basis formula from components in another specified frame. 
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the components are required; if none is provided, 
          the components are assumed to refer to the manifold's default frame
        - ``from_frame`` -- (default: None) name of the frame from which the
          required components are computed, via the tensor change-of-basis 
          formula, if they are not known already in the frame ``frame_name``; 
          if none, a frame name is picked in ``self.components``.
 
        OUTPUT: 
        
        - components in the frame ``frame_name``, as an instance of the 
          class :class:`Components` 
        
        EXAMPLES:
        
        Components of a tensor field of type (1,1) on a 3-dimensional manifold::
    
            sage: m = Manifold(3, 'M')
            sage: c_xyz = Chart(m, 'x y z', 'xyz-coord')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: t = TensorField(m, 1, 1, 'T')
            sage: t.components  # no components defined yet
            {}
            sage: t.set_comp('e') # the first call creates the components
            2-indices components w.r.t. the vector frame 'e'
            sage: t.components
            {'e': 2-indices components w.r.t. the vector frame 'e'}
            
        For the manifold's default frame, the name of the frame is optional::
        
            sage: t.comp()
            2-indices components w.r.t. the vector frame 'e'

        The access to each component can be performed either via the method 
        :meth:`comp` or, in the case of the default frame, via the square
        brackets::
        
            sage: t[0,0] = -2
            sage: t[0,0]
            -2
            sage: t[0,0] is t.comp()[0,0]
            True
            sage: t[0,0] is t.comp('e')[0,0] 
            True
        
        The square brackets return only the components with respect to 
        the default frame::
        
            sage: f = VectorFrame(m, 'f')  # another frame, not the default one
            sage: t.set_comp('f', delete_others=False)[0,0] = -2
            sage: bool( t[0,0] == t.comp('f')[0,0] ) # in the present case, the values are equal :
            True
            sage: t[0,0] is t.comp('f')[0,0] # but not the components:
            False

        Change of components for a vector field on a 2-dimensional manifold::
        
            sage: m = Manifold(2)
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: change_basis = AutomorphismField(m)
            sage: change_basis[:] = [[1,2],[3,4]]
            sage: f = e.new_frame(change_basis, 'f') 
            sage: v = VectorField(m)
            sage: v[0], v[1] = (1,2)   # components set in the default frame (e)
            sage: v.comp('f') # forces the computation of components in frame f
            1-index components w.r.t. the vector frame 'f'
            sage: v.comp('f')[:]
            [0, 1/2]

        Change of components for a type (1,1) tensor field::
        
            sage: t = TensorField(m, 1, 1)
            sage: t[0,0], t[0,1] = (-1, 2)  # components set in the default frame (e)
            sage: t[1,0], t[1,1] = (1, -3)
            sage: t.comp('f')   # forces the computation of components in frame f
            2-indices components w.r.t. the vector frame 'f'
            sage: t_e = matrix([[t[i,j].expr() for j in range(2)] for i in range(2)]) ; t_e
            [-1  2]
            [ 1 -3]
            sage: t_f = matrix([[t.comp('f')[i,j].expr() for j in range(2)] for i in range(2)]) ; t_f
            [ -18  -22]
            [23/2   14]
            
        As a check, one can verify the tensor change-of-basis formula at 
        the matrix level::
        
            sage: p = m.frame_change('e','f')
            sage: p_mat = matrix([[p[i,j].expr() for j in range(2)] for i in range(2)])
            sage: t_f == p_mat.inverse() * t_e * p_mat
            True

        From the Cartesian components to the spherical one for a vector on 
        `\RR^3`::
        
            sage: m = Manifold(3, 'R3', '\RR^3', start_index=1)
            sage: c_cart = Chart(m, 'x y z', 'cart')
            sage: c_spher = Chart(m, r'r:positive th:\theta  ph:\phi','spher')
            sage: assume(th>=0); assume(th<=pi)
            sage: spher_to_cart = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            sage: cart_to_spher = spher_to_cart.set_inverse(sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y, x))
            Check of the inverse coordinate transformation:
              r == r
              th == arctan2(r*sin(th), r*cos(th))
              ph == arctan2(r*sin(ph)*sin(th), r*sin(th)*cos(ph))
              x == x
              y == y
              z == z
            sage: v = VectorField(m)
            sage: v[1], v[2], v[3] = var('vx vy vz') # components set in the Cartesian frame
            sage: v.comp('spher_b') # computation of the components in the spherical frame
            1-index components w.r.t. the coordinate basis 'spher_b' (d/dr,d/dth,d/dph)
            sage: v.comp('spher_b')[1]  # v^r
            (vx*x + vy*y + vz*z)/sqrt(x^2 + y^2 + z^2)
            sage: v.comp('spher_b')[1, 'spher'] # v^r expressed in terms of spherical coord.
            vx*sin(th)*cos(ph) + vy*sin(ph)*sin(th) + vz*cos(th)
            sage: v.comp('spher_b')[2]  # v^\theta
            -sqrt(x^2 + y^2)*(vz*x^2 + vz*y^2 - (vx*x + vy*y)*z)/(x^4 + 2*x^2*y^2 + y^4 + (x^2 + y^2)*z^2)
            sage: v.comp('spher_b')[2, 'spher'] # v^\theta expressed in terms of spherical coord.
            (vx*cos(ph)*cos(th) + vy*sin(ph)*cos(th) - vz*sin(th))/r
            sage: v.comp('spher_b')[3]  # v^\phi
            -(vx*y - vy*x)/(x^2 + y^2)
            sage: v.comp('spher_b')[3, 'spher'] # v^\phi expressed in terms of spherical coord.
            -(vx*sin(ph) - vy*cos(ph))/(r*sin(th))

        The components in the orthonormal basis associated with the spherical 
        coordinates::
        
            sage: change_basis = AutomorphismField(m)
            sage: change_basis.set_comp('spher_b')[1,1, 'spher'] = 1            
            sage: change_basis.set_comp('spher_b')[2,2, 'spher'] = 1/r          
            sage: change_basis.set_comp('spher_b')[3,3, 'spher'] = 1/(r*sin(th))
            sage: e_spher = m.frame('spher_b').new_frame(change_basis, 'e_spher')
            sage: v.comp('e_spher', 'spher_b') # components in 'e_spher' from those in 'spher_b'
            1-index components w.r.t. the vector frame 'e_spher'
            sage: v.comp('e_spher')[1, 'spher']  # v^{(r)}             
            vx*sin(th)*cos(ph) + vy*sin(ph)*sin(th) + vz*cos(th)
            sage: v.comp('e_spher')[2, 'spher']  # v^{(\theta)}
            vx*cos(ph)*cos(th) + vy*sin(ph)*cos(th) - vz*sin(th)
            sage: v.comp('e_spher')[3, 'spher']  # v^{(\phi)}
            -vx*sin(ph) + vy*cos(ph)

        """
        if frame_name is None: 
            frame_name = self.manifold.def_frame.name
        if frame_name not in self.components:
            # the components must be computed from those in the frame from_frame
            manif = self.manifold
            if from_frame is None: 
                from_frame = self.pick_a_frame()
            if from_frame not in self.components:
                raise ValueError("The tensor components are not known in the " +
                                 "frame "+ from_frame)
            (n_con, n_cov) = self.tensor_type
            if n_cov > 0:
                if (from_frame, frame_name) not in manif.frame_changes:
                    raise ValueError("The change-of-basis matrix from the " + 
                                     "frame " + from_frame + " to " + frame_name 
                                     + " has not been set.")
                pp = \
                  manif.frame_changes[(from_frame, frame_name)].comp(from_frame)
                # pp not used if n_cov = 0 (pure contravariant tensor)
            if n_con > 0:
                if (frame_name, from_frame) not in manif.frame_changes:
                    raise ValueError("The change-of-basis matrix from the " + 
                                     "frame " + frame_name + " to " + from_frame 
                                     + " has not been set.")
                ppinv = \
                  manif.frame_changes[(frame_name, from_frame)].comp(from_frame)
                # ppinv not used if n_con = 0 (pure covariant tensor)
            old_comp = self.components[from_frame]
            new_comp = self._new_comp(frame_name)
            si = manif.sindex
            rank = self.rank
            # loop on the new components:
            for ind_new in new_comp.non_redundant_index_generator(): 
                # Summation on the old components multiplied by the proper 
                # change-of-basis matrix elements (tensor formula): 
                res = 0 
                for ind_old in manif.index_generator(rank): 
                    t = old_comp[[ind_old]]
                    for i in range(n_con): # loop on contravariant indices
                        t *= ppinv[[ind_new[i], ind_old[i]]]
                    for i in range(n_con,rank):  # loop on covariant indices
                        t *= pp[[ind_old[i], ind_new[i]]]
                    res += t
                new_comp[ind_new] = res
            self.components[frame_name] = new_comp
            # end of case where the computation was necessary
        return self.components[frame_name]
        
        
    def set_comp(self, frame_name=None, delete_others=True):
        r"""
        Return the components in a given frame for assignment.
        
        NB: the default value of ``delete_others`` ensures that the components 
        with respect to other frames are deleted, in order to avoid any 
        inconsistency. 
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the components are defined; if none is provided, 
          the components are assumed to refer to the manifold's default frame.
        - ``delete_others`` -- (default: True) determines whether the components
          with respect to vector frames different from that specified by 
          ``frame_name`` are deleted or not. 
          
        .. WARNING::
        
            Setting ``delete_others`` to False is at the responsability of the 
            user, who must make sure that the various sets of components are
            consistent between each other. 
         
        OUTPUT: 
        
        - components in the given frame, as an instance of the 
          class :class:`Components`; if such components did not exist
          previously, they are created.  
        
        EXAMPLES:
        
        Components of a tensor field of type (1,1) on a 3-dimensional manifold::
    
            sage: m = Manifold(3, 'M')
            sage: c_xyz = Chart(m, 'x y z', 'xyz-coord')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: t = TensorField(m, 1, 1, 'T')
            sage: t.components  # no components defined yet
            {}
            sage: t.set_comp('e') # the first call creates the components
            2-indices components w.r.t. the vector frame 'e'
            sage: t.components
            {'e': 2-indices components w.r.t. the vector frame 'e'}
            
        For the manifold's default frame, the name of the frame is optional::
        
            sage: t.set_comp()
            2-indices components w.r.t. the vector frame 'e'

        The access to each component can be performed either via the method 
        :meth:`set_comp` or, in the case of the default frame, via the square
        brackets::
        
            sage: t.set_comp()[0,0] = -2
            sage: t[0,0]
            -2
            sage: t[0,0] is t.comp()[0,0]
            True
            sage: t[0,0] is t.comp('e')[0,0] 
            True
        
        The square brackets return only the components with respect to 
        the default frame::
        
            sage: f = VectorFrame(m, 'f')  # another frame, not the default one
            sage: t.set_comp('f', delete_others=False)[0,0] = -2
            sage: bool( t[0,0] == t.comp('f')[0,0] ) # in the present case, the values are equal :
            True
            sage: t[0,0] is t.comp('f')[0,0] # but not the components:
            False

        """
        if frame_name is None: frame_name = self.manifold.def_frame.name
        if frame_name not in self.components:
            if frame_name not in self.manifold.frames:
                raise ValueError("The vector frame " + frame_name +
                                 " has not been defined on the " + 
                                 str(self.manifold))
            self.components[frame_name] = self._new_comp(frame_name)
        self._del_derived() # deletes the derived quantities
        if delete_others:
            self.del_other_comp(frame_name)
        return self.components[frame_name]


    def del_other_comp(self, frame_name=None):
        r"""
        Delete all the components but those corresponding to ``frame_name``.
        
        """
        if frame_name is None: frame_name = self.manifold.def_frame.name
        if frame_name not in self.components:
            raise ValueError("The components w.r.t. the vector frame " + 
                             frame_name + " have not been defined.")
        to_be_deleted = []
        for other_frame_name in self.components:
            if other_frame_name != frame_name:
                to_be_deleted.append(other_frame_name)
        for other_frame_name in to_be_deleted:
            del self.components[other_frame_name]

    def __getitem__(self, indices):
        r"""
        Return the component w.r.t. the manifold default frame corresponding 
        to the given indices. 

        INPUT:
        
        - ``indices`` -- list of indices
    
        """
        return self.comp()[indices]
        
    def __setitem__(self, indices, value):
        r"""
        Set the component w.r.t. the manifold default frame corresponding 
        to the given indices.

        INPUT:
        
        - ``indices`` -- list of indices
    
        """        
        self.set_comp()[indices] = value
        
    def __call__(self, *arg):
        r"""
        The tensor field acting on linear forms and vectors as a multilinear 
        map.
        
        INPUT:
        
        - ``*arg`` -- list of k 1-forms and l vectors, self being a tensor field
          of type (k,l). 
          
        """
        from diffform import OneForm
        from vectorfield import VectorField
        from scalarfield import ScalarField
        # Consistency checks:
        p = len(arg)
        if p != self.rank:
            raise TypeError(str(self.rank) + " arguments must be provided.")
        for i in range(self.tensor_type[0]):
            if not isinstance(arg[i], OneForm):
                raise TypeError("The argument no. " + str(i+1) + 
                                " must be a 1-form.")
        for i in range(self.tensor_type[0],p):
            if not isinstance(arg[i], VectorField):
                raise TypeError("The argument no. " + str(i+1) + 
                                " must be a vector field.")
        manif = self.manifold
        for i in range(p):
            if arg[i].manifold != manif:
                raise ValueError("The argument no. " + str(i+1) + 
                                 "is not defined on the same manifold as the " +
                                 "tensor field.")

        # Does each object have components on the manifold's default frame ? 
        def_frame_name = manif.def_frame.name
        frame_name = None
        if def_frame_name in self.components:
            frame_name = def_frame_name
            for i in range(p):
                if def_frame_name not in arg[i].components:
                    frame_name = None
                    break
        if frame_name is None:
            # Search for a common frame
            for k in self.components:
                frame_name = k
                for i in range(p):
                    if k not in arg[i].components:
                        frame_name = None
                        break
                if frame_name is not None: # common frame found ! 
                    break
            if frame_name is None:
                raise ValueError("No common frame for the components.")
        t = self.components[frame_name]
        v = [arg[i].components[frame_name] for i in range(p)]
        
        res = 0
        for ind in manif.index_generator(self.rank):
            prod = t[[ind]]
            for i in range(p):
                prod *= v[i][[ind[i]]]
            res += prod
        # Name of the output:
        res_name = None
        if self.name is not None:
            res_name = self.name + "("
            for i in range(p-1):
                if arg[i].name is not None:
                    res_name += arg[i].name + ","
                else:
                    res_name = None
                    break
            if res_name is not None:
                if arg[p-1].name is not None:
                    res_name += arg[p-1].name + ")"
                else:
                    res_name = None
        res.name = res_name       
        # LaTeX symbol of the output:
        res_latex = None
        if self.latex_name is not None:
            res_latex = self.latex_name + r"\left("
            for i in range(p-1):
                if arg[i].latex_name is not None:
                    res_latex += arg[i].latex_name + ","
                else:
                    res_latex = None
                    break
            if res_latex is not None:
                if arg[p-1].latex_name is not None:
                    res_latex += arg[p-1].latex_name + r"\right)"
                else:
                    res_latex = None
        res.latex_name = res_latex
        return res

    def common_frame(self, other):
        r"""
        Find a common vector frame for the components of ``self`` and 
        ``other``. 
        
        In case of multiple common frames, the manifold's default frame is 
        privileged. 
        If the current components of ``self`` and ``other`` are all relative to
        different frames, a common frame is searched by performing a component
        transformation, via the transformations listed in 
        ``self.manifold.frame_changes``, still privileging transformations to 
        the manifold's default frame.
        
        INPUT:
        
        - ``other`` -- a tensor field
        
        OUPUT:
        
        - name of the common frame; if a common frame is not found, None is 
          returned. 
        
        """
        # Compatibility checks:
        if not isinstance(other, TensorField):
            raise TypeError("The argument must be a tensor field.")
        manif = self.manifold
        if other.manifold != manif:
            raise TypeError("The two tensor fields are not defined on the " +
                            "same manifold.")
        #
        # NB: in all that follows, the denomination "frame" stands actually for
        # a frame name. 
        def_frame = manif.def_frame.name
        #
        # 1/ Search for a common frame among the existing components, i.e. 
        #    without performing any component transformation. 
        #    -------------------------------------------------------------
        if def_frame in self.components and def_frame in other.components:
            return def_frame # the manifold's default frame is privileged
        for sframe in self.components:
            if sframe in other.components:
                return sframe
        #
        # 2/ Search for a common frame via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least 
        # one component transformation to get a common frame
        if def_frame in self.components:
            for oframe in other.components:
                if (oframe, def_frame) in manif.frame_changes:
                    other.comp(def_frame, from_frame=oframe)
                    return def_frame
        if def_frame in other.components:
            for sframe in self.components:
                if (sframe, def_frame) in manif.frame_changes:
                    self.comp(def_frame, from_frame=sframe)
                    return def_frame
        # If this point is reached, then def_frame cannot be a common frame
        # via a single component transformation
        for sframe in self.components:
            for oframe in other.components:
                if (oframe, sframe) in manif.frame_changes:
                    other.comp(sframe, from_frame=oframe)
                    return sframe
                if (sframe, oframe) in manif.frame_changes:
                    self.comp(oframe, from_frame=sframe)
                    return oframe
        #
        # 3/ Search for a common frame via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at two
        # component transformation to get a common frame
        for sframe in self.components:
            for oframe in other.components:
                if (sframe, def_frame) in manif.frame_changes and \
                   (oframe, def_frame) in manif.frame_changes:
                    self.comp(def_frame, from_frame=sframe)
                    other.comp(def_frame, from_frame=oframe)
                    return def_frame
                for frame in manif.frames:
                    if (sframe, frame) in manif.frame_changes and \
                       (oframe, frame) in manif.frame_changes:
                        self.comp(frame, from_frame=sframe)
                        other.comp(frame, from_frame=oframe)
                        return frame
        #
        # If this point is reached, no common frame could be found, even at 
        # the price of component transformations:
        return None
    
    def pick_a_frame(self):
        r"""
        Return a vector frame in which the tensor components are defined. 
        
        The manifold's default frame is privileged. 

        OUTPUT:
        
        - string containing the name of the vector frame

        """
        if self.manifold.def_frame.name in self.components:
            return self.manifold.def_frame.name  # the default frame is privileged
        else:
            # a frame is picked arbitrarily:
            return self.components.items()[0][0]  

    def pick_a_coord_basis(self):
        r"""
        Return a coordinate basis in which the tensor components are expressed. 
        
        The coordinate basis associated to the manifold's default chart is 
        privileged. 
        If none of the current components is relative to a coordinate basis, 
        a coordinate basis is searched by performing a component transformation, 
        via the transformations listed in ``self.manifold.frame_changes``, still
        privileging transformations to the coordinate basis associated to the 
        manifold's default frame.
        
        OUTPUT:
        
        - string containing the name of the coordinate basis fulfilling the 
          above requirements.

        """
        # 1/ Has the tensor already some components in the default coordinate
        #    basis ?
        manif = self.manifold
        def_coordb = manif.def_chart.name + "_b"
        if def_coordb in self.components:
            return def_coordb
        # 2/ If not, has the tensor already some components in any coordinate
        #    basis ?
        for frame in self.components:
            if frame[-2:] == "_b":
                return frame
        # 3/ If not, let us try to perform a component transformation to the 
        #    default coordinate frame:
        for frame in self.components:
            if (frame, def_coordb) in manif.frame_changes:
                self.comp(def_coordb, from_frame=frame)
                return def_coordb
        # 4/ If unsuccessfull, let us try to perform a component transformation
        #    to any coordinate frame:
        for chart in manif.atlas:
            if chart != manif.def_chart: # the case def_chart is treated in 3/
                coordb = chart + "_b"
                for frame in self.components:
                    if (frame, coordb) in manif.frame_changes:
                        self.comp(coordb, from_frame=frame)
                        return coordb
        # If this point is reached, it has not been possible to find some 
        # components in a coordinate basis, even at the price of a component 
        # transformation:
        return None
        

    def is_zero(self):
        r""" 
        Return True if the tensor field is zero and False otherwise.

        EXAMPLES:
        
        A just-created tensor field with components initialized with 
        :meth:`set_comp` is zero::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: t = TensorField(m, 1, 1, 'T')
            sage: t.set_comp()
            2-indices components w.r.t. the coordinate basis 'coord_xy_b' (d/dx,d/dy)
            sage: t.is_zero()
            True
            sage: t[:]  # indeed:
            [0 0]
            [0 0]
            sage: t[1,1] = 2
            sage: t.is_zero()
            False
            sage: t[1,1] = 0
            sage: t.is_zero()
            True

        It is equivalent to use the operator == to compare to zero::
        
            sage: t == 0
            True
            sage: t != 0
            False
            sage: t[1,1] = 2
            sage: t == 0
            False
            sage: t != 0
            True

        Comparing to a nonzero number is meaningless::
    
            sage: t == 1
            Traceback (most recent call last):
            ...
            TypeError: Cannot compare a tensor field to a number.
        
        """
        frame_name = self.pick_a_frame()
        return self.components[frame_name].is_zero()

    
    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a tensor field or 0
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if self.rank == 0:
            raise NotImplementedError("ScalarField.__eq__ should be called " + 
                                      "instead.")
        if isinstance(other, (int, Integer)): # other should be 0
            if other == 0:
                return self.is_zero()
            else:
                raise TypeError("Cannot compare a tensor field to a number.")
        else: # other is another tensor field
            if not isinstance(other, TensorField):
                raise TypeError("An instance of TensorField is expected.")
            #!# base manifolds should be compared
            if other.tensor_type != self.tensor_type:
                return False
            frame_name = self.common_frame(other)
            if frame_name is None:
                raise ValueError("No common frame for the comparison.")
            return bool(self.components[frame_name] == 
                        other.components[frame_name])

    def __ne__(self, other):
        r"""
        Inequality operator. 
        
        INPUT:
        
        - ``other`` -- a tensor field or 0
        
        OUTPUT:
        
        - True if ``self`` is different from ``other``,  or False otherwise
        
        """
        return not self.__eq__(other)


    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - an exact copy of ``self``
    
        """
        result = self._new_instance()
        for frame_name in self.components:
            result.components[frame_name] = + self.components[frame_name]
        if self.name is not None:
            result.name = '+' + self.name 
        if self.latex_name is not None:
            result.latex_name = '+' + self.latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the tensor field `-T`, where `T` is ``self``
    
        """
        result = self._new_instance()
        for frame_name in self.components:
            result.components[frame_name] = - self.components[frame_name]
        if self.name is not None:
            result.name = '-' + self.name 
        if self.latex_name is not None:
            result.latex_name = '-' + self.latex_name
        return result


    def __add__(self, other):
        r"""
        Tensor addition. 
        
        INPUT:
        
        - ``other`` -- a tensor field, of the same type as ``self``
        
        OUPUT:
        
        - the tensor field resulting from the addition of ``self`` and 
          ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, TensorField):
            raise TypeError("For the addition, other must be a tensor field.")
        if other.tensor_type != self.tensor_type:
            raise TypeError("The two tensor fields are not of the same type.")
        frame_name = self.common_frame(other)
        if frame_name is None:
            raise ValueError("No common frame for the addition.")
        comp_result = self.components[frame_name] + other.components[frame_name]
        result = tensor_field_from_comp(comp_result, self.tensor_type)
        if self.name is not None and other.name is not None:
            result.name = self.name + '+' + other.name
        if self.latex_name is not None and other.latex_name is not None:
            result.latex_name = self.latex_name + '+' + other.latex_name
        return result
    
    def __radd__(self, other):
        r"""
        Addition on the left with ``other``. 
        
        """
        return self.__add__(other)


    def __sub__(self, other):
        r"""
        Tensor subtraction. 
        
        INPUT:
        
        - ``other`` -- a tensor field, of the same type as ``self``
        
        OUTPUT:
        
        - the tensor field resulting from the subtraction of ``other`` from 
          ``self``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, TensorField):
            raise TypeError("For the subtraction, other must be a tensor " + 
                            "field.")
        if other.tensor_type != self.tensor_type:
            raise TypeError("The two tensor fields are not of the same type.")
        frame_name = self.common_frame(other)
        if frame_name is None:
            raise ValueError("No common frame for the subtraction.")
        comp_result = self.components[frame_name] - other.components[frame_name]
        result = tensor_field_from_comp(comp_result, self.tensor_type)
        if self.name is not None and other.name is not None:
            result.name = self.name + '-' + other.name
        if self.latex_name is not None and other.latex_name is not None:
            result.latex_name = self.latex_name + '-' + other.latex_name
        return result

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``. 
        
        """
        return (-self).__add__(other)


    def __mul__(self, other):
        r"""
        Tensor product. 
        """
        from scalarfield import ScalarField
        from utilities import format_mul_txt, format_mul_latex
        if isinstance(other, TensorField):
            if isinstance(other, ScalarField):
                raise NotImplementedError("Product with a scalar field not " + 
                                          "implemented yet.")
            frame_name = self.common_frame(other)
            if frame_name is None:
                raise ValueError("No common frame for the tensor product.")
            comp_prov = self.components[frame_name] * other.components[frame_name]
            # Reordering of the contravariant and covariant indices:
            k1, l1 = self.tensor_type
            k2, l2 = other.tensor_type
            if l1 != 0:
                comp_result = comp_prov.swap_adjacent_indices(k1, self.rank, 
                                                              self.rank+k2)
            else:
                comp_result = comp_prov  # no reordering is necessary
            result = tensor_field_from_comp(comp_result, (k1+k2, l1+l2))
            result.name = format_mul_txt(self.name, '*', other.name)
            result.latex_name = format_mul_latex(self.latex_name, r'\otimes ', 
                                                 other.latex_name)
            return result
        else:
            # multiplication by a scalar: 
            result = self._new_instance()
            for frame_name in self.components:
                result.components[frame_name] = other * self.components[frame_name]
            return result


    def __rmul__(self, other):
        r"""
        Multiplication on the left by ``other``. 
        
        """
        if isinstance(other, TensorField):
            raise NotImplementedError("Left tensor product not implemented.")
        # Left multiplication by a scalar: 
        result = self._new_instance()
        for frame_name in self.components:
            result.components[frame_name] = other * self.components[frame_name]
        return result

    def __div__(self, other):
        r"""
        Division (by a scalar)
        """
        result = self._new_instance()
        for frame_name in self.components:
            result.components[frame_name] = self.components[frame_name] / other
        return result
        

    def self_contract(self, pos1, pos2):
        r""" 
        Contraction on two slots of the tensor field. 
        
        INPUT:
            
        - ``pos1`` -- position of the first index for the contraction
        - ``pos2`` -- position of the second index for the contraction
          
        OUTPUT:
        
        - tensor field resulting from the (pos1, pos2) contraction
       
        EXAMPLES:
        
        Contraction on the two slots of a type (1,1) tensor field::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz-coord')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: a = EndomorphismField(m)  # type (1,1) tensor field
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: a.self_contract(0,1).expr()  # contraction of slot 0 with slot 1
            15
            sage: a.self_contract(1,0).expr()  # the order of the slots does not matter
            15
            
        The contraction on two slots having the same tensor type cannot occur::
        
            sage: b = TensorField(m, 2, 0) # type (2,0) tensor field
            sage: b[:] =  [[1,2,3], [4,5,6], [7,8,9]]
            sage: b.self_contract(0,1)  # both slots are contravariant
            Traceback (most recent call last):
            ...
            IndexError: Contraction on two contravariant indices is not allowed.

        The contraction either perserves or destroys the symmetries::
        
            sage: b = DiffForm(m, 2)
            sage: b[1,2], b[1,3], b[2,3] = (3, 2, 1)
            sage: t = a*b ; t
            tensor field of type (1,3) on the 3-dimensional manifold 'M'
            sage: # by construction, t is a tensor field antisymmetric w.r.t. its last two slots:
            sage: t.symmetries() 
            no symmetry;  antisymmetry: (2, 3)
            sage: s = t.self_contract(0,1) ; s   # contraction on the first two slots
            2-form on the 3-dimensional manifold 'M'
            sage: s.symmetries()    # the antisymmetry is preserved
            no symmetry;  antisymmetry: (0, 1)
            sage: s[:]
            [  0  45  30]
            [-45   0  15]
            [-30 -15   0]
            sage: s == 15*b  # check
            True
            sage: s = t.self_contract(0,2) ; s   # contraction on the first and third slots
            tensor field of type (0,2) on the 3-dimensional manifold 'M'
            sage: s.symmetries()  # the antisymmetry has been destroyed by the above contraction:
            no symmetry;  no antisymmetry
            sage: s[:]  # indeed:
            [-26  -4   6]
            [-31  -2   9]
            [-36   0  12]
            sage: s[:] == matrix( [[sum(t[k,i,k,j].expr() for k in range(1,4)) for j in range(1,4)] for i in range(1,4)] )  # check 
            True
            
        """
        # The indices at pos1 and pos2 must be of different types: 
        k_con = self.tensor_type[0]
        l_cov = self.tensor_type[1]
        if pos1 < k_con and pos2 < k_con:
            raise IndexError("Contraction on two contravariant indices is " +
                             "not allowed.")
        if pos1 >= k_con and pos2 >= k_con:
            raise IndexError("Contraction on two covariant indices is " +
                             "not allowed.")
        # Frame selection for the computation: 
        if self.manifold.def_frame.name in self.components:
            frame_name = self.manifold.def_frame.name
        else: # a frame is picked arbitrarily:
            frame_name = list(self.components)[0]
        cmp_res = self.components[frame_name].self_contract(pos1, pos2)
        if self.rank == 2:
            return cmp_res   #!# scalar case
        else:
            return tensor_field_from_comp(cmp_res, (k_con-1, l_cov-1))
      

    def contract(self, pos1, other, pos2):
        r""" 
        Contraction with another tensor field. 
        
        INPUT:
            
        - ``pos1`` -- position of the first index (in ``self``) for the 
          contraction
        - ``other`` -- the tensor field to contract with
        - ``pos2`` -- position of the second index (in ``other``) for the 
          contraction
          
        OUTPUT:
        
        - tensor field resulting from the (pos1, pos2) contraction
       
        EXAMPLES:
        
        Contraction of a 1-form and a vector field on a 2-dimensional 
        manifold::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: a = OneForm(m)
            sage: a[:] = (1,2)
            sage: v = VectorField(m)
            sage: v[:] = (-4,3)
            sage: a.contract(0, v, 0).expr()   # contraction of a 1-form and a vector
            2
            sage: # same result as the 1-form applied to the vector:
            sage: a.contract(0, v, 0) == a(v)  
            True
            sage: v.contract(0, a, 0).expr()
            2
            sage: # same result as the vector considered as a type (1,0) tensor acting on the 1-form:
            sage: v.contract(0, a, 0) == v(a)
            True

        Contraction with a type (1,1) tensor field::
        
            sage: t = EndomorphismField(m)  # a type (1,1) tensor field
            sage: t[:] = [[1,2], [3,4]]
            sage: s = t.contract(1, v, 0) ; s
            vector field on the 2-dimensional manifold 'M'
            sage: s[:]
            [2, 0]
            sage: # same result as the endormorphism applied to the vector:
            sage: s == t(v)
            True
            sage: s = a.contract(0, t, 0) ; s
            1-form on the 2-dimensional manifold 'M'
            sage: s[:]
            [7, 10]
            sage: a[1]*t[1,1] + a[2]*t[2,1], a[1]*t[1,2] + a[2]*t[2,2] # check
            (7, 10)

        """
        if not isinstance(other, TensorField):
            raise TypeError("For the contraction, other must be a tensor " + 
                            "field.")
        k1, l1 = self.tensor_type
        k2, l2 = other.tensor_type
        if pos1 < k1 and pos2 < k2:
            raise TypeError("Contraction not possible: the two index " + 
                            "positions are both contravariant.")
        if pos1 >= k1 and pos2 >= k2:
            raise TypeError("Contraction not possible: the two index " + 
                            "positions are both covavariant.")
        frame_name = self.common_frame(other)
        if frame_name is None:
            raise ValueError("No common frame for the contraction.")
        cmp_res = self.components[frame_name].contract(pos1, 
                                            other.components[frame_name], pos2)
        # reordering of the indices to have all contravariant indices first:
        if k2 > 1:
            if pos1 < k1:
                cmp_res = cmp_res.swap_adjacent_indices(k1-1, k1+l1-1, k1+l1+k2-1)
            else:
                cmp_res = cmp_res.swap_adjacent_indices(k1, k1+l1-1, k1+l1+k2-2)
        type_res = (k1+k2-1, l1+l2-1)
        if type_res == (0, 0):
            return cmp_res  #!# scalar case
        else:
            return tensor_field_from_comp(cmp_res, type_res)

    def symmetrize(self, pos=None, frame_name=None):
        r"""
        Symmetrization over some arguments.
        
        INPUT:
        
        - ``pos`` -- (default: None) list of argument positions involved in the 
          symmetrization (with the convention position=0 for the first 
          argument); if none, the symmetrization is performed over all the 
          arguments
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame with respect to which the component computation is to be
          performed; if none, the manifold's default frame will be used if
          the tensor field has already components in it, otherwise another frame
          w.r.t. which the tensor field has components will be picked
                  
        OUTPUT:
        
        - the symmetrized tensor field
          
        EXAMPLES:
        
        Symmetrization of a rank-2 tensor on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: t = TensorField(m, 0, 2, 't')                       
            sage: t[:] = [[x, 2*y, 3*z], [4*x,5*y,6*z], [7*x,8*y,9*z]]
            sage: s = t.symmetrize() ; s
            field of symmetric bilinear forms on the 3-dimensional manifold 'M'
            sage: t[:], s[:]
            (
            [  x 2*y 3*z]  [            x       2*x + y 7/2*x + 3/2*z]
            [4*x 5*y 6*z]  [      2*x + y           5*y     4*y + 3*z]
            [7*x 8*y 9*z], [7/2*x + 3/2*z     4*y + 3*z           9*z]
            )
            sage: s == t.symmetrize([0,1])
            True
        
        Symmetrization is meaningfull only on arguments of the same type::
        
            sage: h = TensorField(m,1,1)  # a type (1,1) tensor field
            sage: h[1,3] = x
            sage: s = h.symmetrize()
            Traceback (most recent call last):
            ...
            TypeError: Symmetrization is meaningfull only on tensor arguments of the same type.
         
        Symmetrization of a rank-3 tensor::
        
            sage: v = VectorField(m)
            sage: v[:] = (-1, 2, -3)
            sage: w = VectorField(m)
            sage: w[:] = (x,y,z)
            sage: u = VectorField(m)
            sage: u[:] = (cos(z), sin(y), cos(x))
            sage: t = u*v*w ; t  # rank-3 tensor constructed via tensor products
            tensor field of type (3,0) on the 3-dimensional manifold 'M'
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: s = t.symmetrize() ; s
            tensor field of type (3,0) on the 3-dimensional manifold 'M'
            sage: s.symmetries()
            symmetry: (0, 1, 2);  no antisymmetry
            sage: s[1,2,3]
            -1/6*(3*y - 2*z)*cos(z) - 1/2*x*sin(y) + 1/3*x*cos(x) - 1/6*y*cos(x) - 1/6*z*sin(y)
            sage: s[1,2,3] == (t[1,2,3]+t[1,3,2]+t[2,3,1]+t[2,1,3]+t[3,1,2]+t[3,2,1])/6  # check
            True
            
        Partial symmetrization::
        
            sage: s = t.symmetrize([0,1]) ; s  # symmetrization over the first two arguments
            tensor field of type (3,0) on the 3-dimensional manifold 'M'
            sage: s.symmetries() 
            symmetry: (0, 1);  no antisymmetry
            sage: s[1,2,3] == (t[1,2,3]+t[2,1,3])/2  # check
            True
            sage: s = t.symmetrize([1,2]) ; s  # symmetrization over the last two arguments   
            tensor field of type (3,0) on the 3-dimensional manifold 'M'
            sage: s.symmetries()                                                       
            symmetry: (1, 2);  no antisymmetry
            sage: s[1,2,3] == (t[1,2,3]+t[1,3,2])/2  # check
            True
            sage: s = t.symmetrize([0,2]) ; s  # the symmetry positions need not to be adjacent
            tensor field of type (3,0) on the 3-dimensional manifold 'M'
            sage: s.symmetries()
            symmetry: (0, 2);  no antisymmetry
            sage: s[1,2,3] == (t[1,2,3]+t[3,2,1])/2  # check
            True
        
        """
        if pos is None:
            pos = range(self.rank)
        # check whether the symmetrization is possible:
        pos_cov = self.tensor_type[0]   # first covariant position 
        pos0 = pos[0]
        if pos0 < pos_cov:  # pos0 is a contravariant position
            for k in range(1,len(pos)):
                if pos[k] >= pos_cov:
                    print "pos[0] is a contravariant position, while pos[" + \
                           str(k)+"] is a covariant position."
                    raise TypeError("Symmetrization is meaningfull only " +
                                    "on tensor arguments of the same type.")
        else:  # pos0 is a covariant position
            for k in range(1,len(pos)):
                if pos[k] < pos_cov:
                    print "pos[0] is a covariant position, while pos[" + \
                           str(k)+"] is a contravariant position."
                    raise TypeError("Symmetrization is meaningfull only " +
                                    "on tensor arguments of the same type.")                
        if frame_name is None:
            frame_name = self.pick_a_frame()
        res_comp = self.components[frame_name].symmetrize(pos)
        return tensor_field_from_comp(res_comp, self.tensor_type)
                

    def antisymmetrize(self, pos=None, frame_name=None):
        r"""
        Antisymmetrization over some arguments.
        
        INPUT:
        
        - ``pos`` -- (default: None) list of argument positions involved in the 
          antisymmetrization (with the convention position=0 for the first 
          argument); if none, the antisymmetrization is performed over all the 
          arguments
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame with respect to which the component computation is to be
          performed; if none, the manifold's default frame will be used if
          the tensor field has already components in it, otherwise another frame
          w.r.t. which the tensor field has components will be picked
                  
        OUTPUT:
        
        - the antisymmetrized tensor field
          
        EXAMPLES:
        
        Antisymmetrization of a rank-2 tensor on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: t = TensorField(m, 0, 2, 't')                       
            sage: t[:] = [[x, 2*y, 3*z], [4*x,5*y,6*z], [7*x,8*y,9*z]]        
            sage: s = t.antisymmetrize() ; s     
            2-form on the 3-dimensional manifold 'M'
            sage: t[:], s[:]
            (
            [  x 2*y 3*z]  [             0       -2*x + y -7/2*x + 3/2*z]
            [4*x 5*y 6*z]  [       2*x - y              0     -4*y + 3*z]
            [7*x 8*y 9*z], [ 7/2*x - 3/2*z      4*y - 3*z              0]
            )
            sage: s == t.antisymmetrize([0,1])
            True

        Antisymmetrization is meaningfull only on arguments of the same type::
        
            sage: h = TensorField(m,1,1)  # a type (1,1) tensor field
            sage: h[1,3] = x
            sage: s = h.antisymmetrize()
            Traceback (most recent call last):
            ...
            TypeError: Antisymmetrization is meaningfull only on tensor arguments of the same type.
         

        Antisymmetrization of a rank-3 tensor::
        
            sage: w = OneForm(m)
            sage: w[:] = (z,y,x) 
            sage: t = w*t ; t  # rank-3 tensor constructed via a tensor product
            tensor field of type (0,3) on the 3-dimensional manifold 'M'
            sage: t[:]
            [[[x*z, 2*y*z, 3*z^2], [4*x*z, 5*y*z, 6*z^2], [7*x*z, 8*y*z, 9*z^2]],
             [[x*y, 2*y^2, 3*y*z], [4*x*y, 5*y^2, 6*y*z], [7*x*y, 8*y^2, 9*y*z]],
             [[x^2, 2*x*y, 3*x*z], [4*x^2, 5*x*y, 6*x*z], [7*x^2, 8*x*y, 9*x*z]]]
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: s = t.antisymmetrize() ; s  # full antisymmetrization: the result is a 3-form:
            3-form on the 3-dimensional manifold 'M'
            sage: s[1,2,3]
            -2/3*x^2 + 3/2*x*y - 11/6*y*z + z^2
            sage: s[1,2,3] == (t[1,2,3]-t[1,3,2]+t[2,3,1]-t[2,1,3]+t[3,1,2]-t[3,2,1])/6  # check
            True
            
        Partial antisymmetrizations of a rank-3 tensor::
            
            sage: s = t.antisymmetrize([0,1]) ; s  # symmetrization over the first two arguments
            tensor field of type (0,3) on the 3-dimensional manifold 'M'
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: s[1,2,3]
            -3/2*y*z + 3*z^2
            sage: s[1,2,3] == (t[1,2,3]-t[2,1,3])/2  # check
            True
            sage: s = t.antisymmetrize([1,2]) ; s  # symmetrization over the last two arguments
            tensor field of type (0,3) on the 3-dimensional manifold 'M'
            sage: s.symmetries()                   
            no symmetry;  antisymmetry: (1, 2)
            sage: s[1,2,3]
            -4*y*z + 3*z^2
            sage: s[1,2,3] == (t[1,2,3]-t[1,3,2])/2  # check
            True
            sage: s = t.antisymmetrize([0,2]) ; s  # symmetrization over the first and last arguments
            tensor field of type (0,3) on the 3-dimensional manifold 'M'
            sage: s.symmetries()                   
            no symmetry;  antisymmetry: (0, 2)
            sage: s[1,2,3]
            -2*x^2 + 3*z^2
            sage: s[1,2,3] == (t[1,2,3]-t[3,2,1])/2  # check
            True
        
        Recovering the exterior product of p-forms by antisymmetrization of the
        tensor product::
        
            sage: a = OneForm(m)
            sage: a[:] = (x*y*z, -z*x, y*z)
            sage: b = OneForm(m)           
            sage: b[:] = (cos(z), sin(x), cos(y))
            sage: a.wedge(b) == 2*(a*b).antisymmetrize()
            True
            sage: c = a.exterior_der() ; c
            2-form on the 3-dimensional manifold 'M'
            sage: b.wedge(c) == 6/(1*2)*(b*c).antisymmetrize()
            True

     
        """
        if pos is None:
            pos = range(self.rank)
        # check whether the antisymmetrization is possible:
        pos_cov = self.tensor_type[0]   # first covariant position 
        pos0 = pos[0]
        if pos0 < pos_cov:  # pos0 is a contravariant position
            for k in range(1,len(pos)):
                if pos[k] >= pos_cov:
                    print "pos[0] is a contravariant position, while pos[" + \
                           str(k)+"] is a covariant position."
                    raise TypeError("Antisymmetrization is meaningfull only " +
                                    "on tensor arguments of the same type.")
        else:  # pos0 is a covariant position
            for k in range(1,len(pos)):
                if pos[k] < pos_cov:
                    print "pos[0] is a covariant position, while pos[" + \
                           str(k)+"] is a contravariant position."
                    raise TypeError("Antisymmetrization is meaningfull only " +
                                    "on tensor arguments of the same type.")                
        if frame_name is None:
            frame_name = self.pick_a_frame()
        res_comp = self.components[frame_name].antisymmetrize(pos)
        return tensor_field_from_comp(res_comp, self.tensor_type)


    def up(self, metric, pos=None):
        r"""
        Compute a metric dual by raising some index with a given metric.
        
        If ``self`` is a tensor field `T` of type `(k,\ell)` and `p` is the 
        position of a covariant index (i.e. `k\leq p < k+\ell`), 
        the output with ``pos`` `=p` is the tensor field `T^\sharp` of type 
        `(k+1,\ell-1)` whose components are

        .. MATH::

            (T^\sharp)^{a_1\ldots a_{k+1}}_{\qquad\quad b_1 \ldots b_{\ell-1}}
            = g^{a_{k+1} i} \, 
            T^{a_1\ldots a_k}_{\qquad\ \  b_1 \ldots b_{p-k} \, i \, b_{p-k+1} \ldots b_{\ell-1}},
            
        `g^{ab}` being the components of the inverse metric. 

        The reverse operation is :meth:`TensorField.down`

        INPUT:
        
        - ``metric`` -- metric `g`, as an instance of :class:`Metric`
        - ``pos`` -- (default: None) position of the index (with the
          convention ``pos=0`` for the first index); if none, the raising is performed
          over all the covariant indices, starting from the first one
         
        OUTPUT:
        
        - the tensor field `T^\sharp` resulting from the index raising operation

        EXAMPLES:
        
        Raising the index of a 1-form results in a vector field::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: g = Metric(m, 'g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-y
            sage: w = OneForm(m)
            sage: w[:] = [-1, 2]
            sage: v = w.up(g) ; v
            vector field on the 2-dimensional manifold 'M'
            sage: v.show()
            ((2*x - 1)*y + 1)/(x^2*y^2 + (x + 1)*y - x - 1) d/dx - (x*y + 2*x + 2)/(x^2*y^2 + (x + 1)*y - x - 1) d/dy
            sage: g.inverse()[:]
            [ (y - 1)/(x^2*y^2 + (x + 1)*y - x - 1)      x*y/(x^2*y^2 + (x + 1)*y - x - 1)]
            [     x*y/(x^2*y^2 + (x + 1)*y - x - 1) -(x + 1)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: w1 = v.down(g) ; w1   # the reverse operation
            1-form on the 2-dimensional manifold 'M'
            sage: w1.show()
            -dx + 2 dy
            sage: w1 == w
            True

        Raising the indices of a tensor field of type (0,2)::

            sage: t = TensorField(m, 0, 2)
            sage: t[:] = [[1,2], [3,4]]
            sage: tu0 = t.up(g, 0) ; tu0  # raising the first index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: tu0[:]
            [  ((3*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1) 2*((2*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1)]
            [    (x*y - 3*x - 3)/(x^2*y^2 + (x + 1)*y - x - 1)   2*(x*y - 2*x - 2)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: tuu0 = tu0.up(g) ; tuu0 # the two indices have been raised, starting from the first one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: tu1 = t.up(g, 1) ; tu1 # raising the second index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: tu1[:]
            [((2*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1) ((4*x + 3)*y - 3)/(x^2*y^2 + (x + 1)*y - x - 1)]
            [  (x*y - 2*x - 2)/(x^2*y^2 + (x + 1)*y - x - 1) (3*x*y - 4*x - 4)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: tuu1 = tu1.up(g) ; tuu1 # the two indices have been raised, starting from the second one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: tuu0 == tuu1 # the order of index raising is important
            False
            sage: tuu = t.up(g) ; tuu # both indices are raised, starting from the first one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: tuu0 == tuu # the same order for index raising has been applied
            True
            sage: tuu1 == tuu # to get tuu1, indices have been raised from the last one, contrary to tuu 
            False
            sage: d0tuu = tuu.down(g, 0) ; d0tuu # the first index is lowered again
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: dd0tuu = d0tuu.down(g) ; dd0tuu  # the second index is then lowered
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: d1tuu = tuu.down(g, 1) ; d1tuu # lowering operation, starting from the last index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: dd1tuu = d1tuu.down(g) ; dd1tuu
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: ddtuu = tuu.down(g) ; ddtuu # both indices are lowered, starting from the last one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: ddtuu == t # should be true
            True
            sage: dd0tuu == t # not true, because of the order of index lowering to get dd0tuu
            False
            sage: dd1tuu == t # should be true
            True

        """
        n_con = self.tensor_type[0] # number of contravariant indices = k 
        if pos is None:
            result = self
            for p in range(n_con, self.rank):
                k = result.tensor_type[0]
                result = result.up(metric, k)
            return result
        if metric.manifold != self.manifold:
            raise TypeError("The metric is not defined on the same manifold " +
                            "as the tensor.")
        if not isinstance(pos, (int, Integer)):
            raise TypeError("The argument 'pos' must be an integer.")
        if pos<n_con or pos>self.rank-1:
            print "pos = ", pos
            raise ValueError("Position out of range.")
        return self.contract(pos, metric.inverse(), 0)
        
        
    def down(self, metric, pos=None):
        r"""
        Compute a metric dual by lowering some index with a given metric.
        
        If ``self`` is a tensor field `T` of type `(k,\ell)` and `p` is the 
        position of a contravariant index (i.e. `0\leq p < k`), the output with
        ``pos`` `=p` is the tensor field `T^\flat` of type `(k-1,\ell+1)` whose 
        components are

        .. MATH::

            (T^\flat)^{a_1\ldots a_{k-1}}_{\qquad\quad b_1 \ldots b_{\ell+1}}
            = g_{b_1 i} \, 
            T^{a_1\ldots a_{p} \, i \, a_{p+1}\ldots a_{k-1}}_{\qquad\qquad\qquad\quad b_2 \ldots b_{\ell+1}},
            
        `g_{ab}` being the components of the metric tensor. 

        The reverse operation is :meth:`TensorField.up`
        
        INPUT:
        
        - ``metric`` -- metric `g`, as an instance of :class:`Metric`
        - ``pos`` -- (default: None) position of the index (with the 
          convention ``pos=0`` for the first index); if none, the lowering is 
          performed over all the contravariant indices, starting from the last one
         
        OUTPUT:
        
        - the tensor field `T^\flat` resulting from the index lowering operation
        
        EXAMPLES:
        
        Lowering the index of a vector field results in a 1-form::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: g = Metric(m, 'g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-y
            sage: v = VectorField(m)
            sage: v[:] = [-1,2]
            sage: w = v.down(g) ; w
            1-form on the 2-dimensional manifold 'M'
            sage: w.show()
            (2*x*y - x - 1) dx + (-(x + 2)*y + 2) dy
            sage: v1 = w.up(g) ; v1  # the reverse operation
            vector field on the 2-dimensional manifold 'M'
            sage: v1 == v
            True

        Lowering the indices of a tensor field of type (2,0)::
        
            sage: t = TensorField(m, 2, 0)
            sage: t[:] = [[1,2], [3,4]]
            sage: td0 = t.down(g, 0) ; td0  # lowering the first index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: td0[:]
            [  3*x*y + x + 1   (x - 3)*y + 3]
            [4*x*y + 2*x + 2 2*(x - 2)*y + 4]
            sage: tdd0 = td0.down(g) ; tdd0 # the two indices have been lowered, starting from the first one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: tdd0[:]
            [      4*x^2*y^2 + 5*(x^2 + x)*y + x^2 + 2*x + 1 2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3]
            [(3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: td1 = t.down(g, 1) ; td1  # lowering the second index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: td1[:]
            [  2*x*y + x + 1   (x - 2)*y + 2]
            [4*x*y + 3*x + 3 (3*x - 4)*y + 4]
            sage: tdd1 = td1.down(g) ; tdd1 # the two indices have been lowered, starting from the second one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: tdd1[:]
            [      4*x^2*y^2 + 5*(x^2 + x)*y + x^2 + 2*x + 1 (3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2]
            [2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: tdd1 == tdd0   # the order of index lowering is important
            False
            sage: tdd = t.down(g) ; tdd  # both indices are lowered, starting from the last one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: tdd[:]
            [      4*x^2*y^2 + 5*(x^2 + x)*y + x^2 + 2*x + 1 (3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2]
            [2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: tdd0 == tdd  # to get tdd0, indices have been lowered from the first one, contrary to tdd 
            False
            sage: tdd1 == tdd  # the same order for index lowering has been applied
            True
            sage: u0tdd = tdd.up(g, 0) ; u0tdd # the first index is raised again
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: uu0tdd = u0tdd.up(g) ; uu0tdd # the second index is then raised
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: u1tdd = tdd.up(g, 1) ; u1tdd  # raising operation, starting from the last index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: uu1tdd = u1tdd.up(g) ; uu1tdd
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: uutdd = tdd.up(g) ; uutdd  # both indices are raised, starting from the first one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: uutdd == t  # should be true
            True
            sage: uu0tdd == t # should be true
            True
            sage: uu1tdd == t # not true, because of the order of index raising to get uu1tdd
            False
 
        """
        n_con = self.tensor_type[0] # number of contravariant indices = k 
        if pos is None:
            result = self
            for p in range(0, n_con):
                k = result.tensor_type[0]
                result = result.down(metric, k-1)
            return result
        if metric.manifold != self.manifold:
            raise TypeError("The metric is not defined on the same manifold " +
                            "as the tensor.")
        if not isinstance(pos, (int, Integer)):
            raise TypeError("The argument 'pos' must be an integer.")
        if pos<0 or pos>=n_con:
            print "pos = ", pos
            raise ValueError("Position out of range.")
        return metric.contract(1, self, pos)

        
    def lie_der(self, vector):
        r"""
        Computes the Lie derivative with respect to a vector field.
        
        The Lie derivative is stored in the dictionary 
        :attr:`_lie_derivatives`, so that there is no need to 
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile. 
        
        INPUT:
        
        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken
          
        OUTPUT:
        
        - the tensor field that is the Lie derivative of ``self`` with respect 
          to ``vector``
        
        EXAMPLES:
        
        Lie derivative of a vector::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: v = VectorField(m, 'v')
            sage: v[:] = (-y, x)
            sage: w = VectorField(m)
            sage: w[:] = (2*x+y, x*y)
            sage: w.lie_der(v)
            vector field on the 2-dimensional manifold 'M'
            sage: w.lie_der(v).show()
            ((x - 2)*y + x) d/dx + (x^2 - y^2 - 2*x - y) d/dy

        The Lie derivative is antisymmetric::
        
            sage: w.lie_der(v) == -v.lie_der(w)
            True
            
        For vectors, it coincides with the commutator::

            sage: f = ScalarField(m, x^3 + x*y^2)
            sage: w.lie_der(v)(f).show()
            (x, y) |--> -(x + 2)*y^3 + 3*x^3 - x*y^2 + 5*(x^3 - 2*x^2)*y
            sage: w.lie_der(v)(f) == v(w(f)) - w(v(f))  # rhs = commutator [v,w] acting on f
            True
            
        Lie derivative of a 1-form::
        
            sage: om = OneForm(m)
            sage: om[:] = (y^2*sin(x), x^3*cos(y))
            sage: om.lie_der(v)
            1-form on the 2-dimensional manifold 'M'
            sage: om.lie_der(v).show()
            (x^3*cos(y) - y^3*cos(x) + 2*x*y*sin(x)) dx + (-x^4*sin(y) - 3*x^2*y*cos(y) - y^2*sin(x)) dy
            
        Check of Cartan identity::
        
            sage: om.lie_der(v) == v.contract(0,om.exterior_der(),0) + (om(v)).exterior_der()
            True
        
        """
        from scalarfield import scalar_field
        from vectorfield import VectorField
        from vectorframe import CoordBasis
        if not isinstance(vector, VectorField):
            raise TypeError("The argument must be a vector field.")
        manif = self.manifold
        if vector.manifold != manif:
            raise TypeError("The vector field and the tensor field are not " + 
                            "defined on the same manifold")
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            #
            # 1/ Search for a common coordinate basis:
            chart_name = self.manifold.def_chart.name
            coord_basis_name = chart_name + '_b'
            if (coord_basis_name not in self.components) or \
               (coord_basis_name not in vector.components):
                chart_name = None
                coord_basis_name = None
                for frame_name in self.components:
                    if (frame_name in vector.components) and \
                       isinstance(manif.frames[frame_name], CoordBasis):
                        coord_basis_name = frame_name
                        chart_name = coord_basis_name[:-2]
                        break
                if coord_basis_name is None:
                    raise TypeError("No common coordinate basis found.")
            #
            # 2/ Component computation:
            tc = self.components[coord_basis_name]
            vc = vector.components[coord_basis_name]
            # the result has the same tensor type and same symmetries as self:
            resc = self._new_comp(coord_basis_name) 
            si = manif.sindex
            nsi = manif.dim + si
            n_con = self.tensor_type[0]
            for ind in resc.non_redundant_index_generator():
                rsum = 0
                for i in range(si, nsi):
                    rsum += vc[[i]].function_chart(chart_name) * \
                           tc[[ind]].function_chart(chart_name).diff(i)
                # loop on contravariant indices:
                for k in range(n_con): 
                    for i in range(si, nsi):
                        indk = list(ind)
                        indk[k] = i  
                        rsum -= tc[[indk]].function_chart(chart_name) * \
                                vc[[ind[k]]].function_chart(chart_name).diff(i)
                # loop on covariant indices:
                for k in range(n_con, self.rank): 
                    for i in range(si, nsi):
                        indk = list(ind)
                        indk[k] = i  
                        rsum += tc[[indk]].function_chart(chart_name) * \
                                vc[[i]].function_chart(chart_name).diff(ind[k])
                resc[[ind]] = scalar_field(rsum)
            #
            # 3/ Final result (the tensor)
            res = self._new_instance()
            res.components[coord_basis_name] = resc
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]
        
#******************************************************************************
    
def tensor_field_from_comp(comp, tensor_type):
    r"""
    Constructs a tensor field from components in a given frame.
    
    The symmetries/antisymmetries on the components' indices are forwarded as
    symmetries/antisymmetries on the tensor's arguments. 
    Moreover, instances of subclasses of :class:`TensorField` are returned when
    appropriate; for instance, if ``tensor_type`` is (1,0), an instance of
    :class:`VectorField` is returned. 
    
    INPUT:
    
    - ``comp`` : a set of components, as an instance of :class:`Components`
    - ``tensor_type`` : the pair `(k,\ell)` defining the tensor type
    
    OUTPUT:
    
    - an instance of :class:`TensorField` or one of its derived classes. 
    
    EXAMPLES:
    
    A symmetric bilinear form constructed from symmetric 2-indices components::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = Chart(m, 'x y z', 'xyz-coord')
        sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
        sage: c = CompFullySym(m, 2)
        sage: c[0,0], c[0,1], c[0,2] = (1,2,3)
        sage: t = tensor_field_from_comp(c, (0,2)) ; t 
        field of symmetric bilinear forms on the 3-dimensional manifold 'M'
        sage: t.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: t[:]
        [1 2 3]
        [2 0 0]
        [3 0 0]
        sage: (t[:] - c[:]).is_zero()
        True
    
    """
    from vectorfield import VectorField
    from rank2field import SymBilinFormField, EndomorphismField
    from diffform import DiffForm, OneForm
    if not isinstance(comp, Components):
        raise TypeError("The first argument must be a set of components.")
    k = tensor_type[0]
    l = tensor_type[1]
    rank = k+l
    if comp.nid != rank:
        raise TypeError("The number of indices is not equal to the tensor " + 
                        "rank.")
    if (k,l) == (1,0):
        result = VectorField(comp.manifold)
    elif (k,l) == (0,1):
        result = OneForm(comp.manifold)
    elif (k,l) == (0,2) and isinstance(comp, CompFullySym):
        result = SymBilinFormField(comp.manifold)
    elif (k,l) == (1,1):
        result = EndomorphismField(comp.manifold)
    elif k == 0 and isinstance(comp, CompFullyAntiSym):
        result = DiffForm(comp.manifold, rank)
    elif isinstance(comp, CompWithSym):
        result = TensorField(comp.manifold, k, l, sym=comp.sym,
                             antisym=comp.antisym)
    else:
        result = TensorField(comp.manifold, k, l)
       
    result.components[comp.frame.name] = comp
    return result
   
