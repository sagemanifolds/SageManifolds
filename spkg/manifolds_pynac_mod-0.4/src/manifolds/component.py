r"""
Components

The class :class:`Components` takes in charge the storage of the components of 
a geometrical entity with respect to a given vector frame. This concerns of 
course tensor components (see :meth:`TensorField.comp`), but also 
non-tensorial quantities, like connection coefficients (see 
:meth:`AffConnection.coef`) or structure coefficients of a vector frame (see
:meth:`VectorFrame.structure_coef`). 

Various subclasses of the class :class:`Components` are

* :class:`CompWithSym` for storing components with symmetries (symmetric and/or 
  antisymmetric indices)
* :class:`CompFullySym` for storing fully symmetric components.
* :class:`CompFullyAntiSym` for storing fully antisymmetric components.
* :class:`KroneckerDelta` for the Kronecker delta symbol. 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES:

    2-indices set of components on a 3-dimensional manifold::

        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = m.chart('x y z')
        sage: m.default_frame()  # the frame associated with the chart (M, (x,y,z))
        coordinate frame (M, (d/dx,d/dy,d/dz))
        sage: c = Components(m.default_frame(), 2) ; c
        2-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))

    The individual components are accessed by providing their indices inside
    double square brackets::

        sage: c[[1,2]] = x+1

    Each individual component is an instance of :class:`ScalarField`::
    
        sage: c[[1,2]]
        scalar field on the 3-dimensional manifold 'M'
        sage: c[[1,2]].expr()
        x + 1

    A short-cut of the above is the use of single square brackets::
    
        sage: c[1,2]
        x + 1
        sage: c[1,2] == c[[1,2]].expr()
        True
        sage: c[1,2] is c[[1,2]].function_chart()
        True
        sage: c[1,2].expr() is c[[1,2]].expr() # the symbolic expression
        True

    In other words, the single square brackets return an instance of 
    :class:`FunctionChart` that is the coordinate function expressing the 
    component in some chart (by default, the manifold's default chart)::
    
        sage: print type(c[1,2])    # single bracket --> FunctionChart
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: print type(c[[1,2]])  # double bracket --> ScalarField
        <class 'sage.geometry.manifolds.scalarfield.ScalarField'>
    
    Expressions in a chart different from the manifold's default one are 
    obtained by specifying the chart as the last argument inside the
    single square brackets::
    
        sage: c_uvw.<u,v,w> = m.chart('u v w')
        sage: xyz_to_uvw = CoordChange(c_xyz, c_uvw, x+y, x-y, x+y+z)
        sage: uvw_to_xyz = xyz_to_uvw.inverse()
        sage: c[1,2, c_uvw]
        1/2*u + 1/2*v + 1
        sage: c[1,2, c_uvw] == c[[1,2]].expr(c_uvw)
        True

    Components that have not been set are initialized to zero::
    
        sage: c[2,1]
        0
        sage: c[[2,1]]
        zero scalar field on the 3-dimensional manifold 'M'

    The list of all components is provided by [:]::
    
        sage: c[:]
        [    0 x + 1     0]
        [    0     0     0]
        [    0     0     0]
    
    Internally, the components are stored as a dictionary (:attr:`_comp`) whose
    keys are the indices; only the non-zero components are stored::

        sage: c._comp
        {(1, 2): scalar field on the 3-dimensional manifold 'M'}

    In case of symmetries, only the non-redundant components are stored::
    
        sage: c = CompFullyAntiSym(m.default_frame(), 2) # an antisymmetric set of 2-indices components
        sage: c[[1,2]] = x+1
        sage: c[2,1]  # c_{21} is not zero...
        -x - 1
        sage: c._comp  # ...by only c_{12} is stored, since c_{21} is deductible from it by antisymmetry:
        {(1, 2): scalar field on the 3-dimensional manifold 'M'}
        sage: c[:]
        [     0  x + 1      0]
        [-x - 1      0      0]
        [     0      0      0]

    
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
from vectorframe import VectorFrame


class Components(SageObject):
    r"""
    Class for storing components with respect to a given vector frame on a 
    differentiable manifold over `\RR`. 
    
    The stored quantities can be tensor components or non-tensorial quantities, 
    such as connection coefficients. The symmetries over some indices are
    treated by subclasses of the class :class:`Components`.
     
    INPUT:
    
    - ``frame`` -- vector frame with respect to which the components are 
      defined
    - ``nb_indices`` -- number of indices 
    
    EXAMPLES:

    3-indices components on a 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy.<x,y> = m.chart('x y') 
        sage: c = Components(m.default_frame(), 3) ; c
        3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy))
        
    When created, components are initialized to zero::
    
        sage: c[:]  # returns the list of all components
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
        sage: c.is_zero() 
        True
        sage: c == 0  # equivalent to c.is_zero() 
        True

    Components with respect to a frame that is not the manifold's default one::
    
        sage: e = VectorFrame(m, 'e')  # a second frame on M
        sage: m.default_frame()
        coordinate frame (M, (d/dx,d/dy))
        sage: c1 = Components(e, 3) ; c1
        3-indices components w.r.t. the vector frame (M, (e_0,e_1))
 
    The basic attributes are :attr:`manifold`, :attr:`domain`, :attr:`frame` and :attr:`nid` 
    (number of indices)::
        
        sage: c.manifold
        2-dimensional manifold 'M'
        sage: c.domain
        2-dimensional manifold 'M'
        sage: c.frame
        coordinate frame (M, (d/dx,d/dy))
        sage: c.nid
        3

    The access to the components is performed by providing the list of indices
    inside square brackets::
    
        sage: c[0,0,0] = 1 ; c[0,0,1] = pi
        sage: c[0,0,0]
        1
        sage: c == 0  # the set of components is no longer zero:
        False
    
    The access to all the components is provided by the operator [:]::
    
        sage: c[:]  # read access
        [[[1, pi], [0, 0]], [[0, 0], [0, 0]]]
        sage: c[:] = [[[2, 4], [3, 6]], [[4, 8], [5, 10]]]  # write access
        sage: c[0,0,1]
        4

    The components can be symbolic expressions::
    
        sage: c[0,0,0] = x + 2*y
        sage: c[0,0,0] / c[0,0,1]
        1/4*x + 1/2*y
        
    The index range depends on the manifold::
    
        sage: m = Manifold(2, 'M', start_index=1)
        sage: c_xy = m.chart('x y') 
        sage: c = Components(m.default_frame(), 3)
        sage: c[0,0,0] = 1
        Traceback (most recent call last):
        ...
        IndexError: Index out of range: 0 not in [1,2]
        sage: c[1,2,1] = 1  # OK 

    Tensor products of set of components are performed by means of the 
    operator \*::
        
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z')
        sage: e = VectorFrame(m, 'e')
        sage: a = Components(e, 1)
        sage: a[0], a[1], a[2] = (0, 1, 2)
        sage: b = Components(e, 1)
        sage: b[0], b[1], b[2] = (3, 4, 5)
        sage: c = a*b ; c
        2-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
        sage: c[:]
        [ 0  0  0]
        [ 3  4  5]
        [ 6  8 10]
        sage: d = a*c ; d
        3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
        sage: d[:]
        [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 0, 0], [3, 4, 5], [6, 8, 10]], [[0, 0, 0], [6, 8, 10], [12, 16, 20]]]

    The tensor product propagates the symmetries::
        
        sage: f = CompFullySym(e,2)
        sage: f[0,0], f[0,1], f[0,2] = (3, 4, 5)
        sage: f[1,1], f[1,2] = (6,7)
        sage: f[2,2] = 8
        sage: c = a*f ; c
        3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2)), with symmetry on the index positions (1, 2)
        sage: a[:]
        [0, 1, 2]
        sage: f[:]
        [3 4 5]
        [4 6 7]
        [5 7 8]
        sage: c[:]
        [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [[3, 4, 5], [4, 6, 7], [5, 7, 8]], [[6, 8, 10], [8, 12, 14], [10, 14, 16]]]
        sage: type(c)
        <class 'sage.geometry.manifolds.component.CompWithSym'>

    """
    def __init__(self, frame, nb_indices):
        if not isinstance(frame, VectorFrame):
            raise TypeError("The first argument must be a vector frame.")
        if not isinstance(nb_indices, (int, Integer)):
            raise TypeError("The second argument must be an integer.")
        if nb_indices < 1:
            raise TypeError("nb_indices must be equal to or larger than 1.")
        self.frame = frame
        self.manifold = frame.manifold
        self.domain = frame.domain
        self.nid = nb_indices
        self._comp = {} # the dictionary of components, with the indices as keys
        
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = str(self.nid)
        if self.nid == 1:
            description += "-index"
        else:
            description += "-indices"
        description += " components w.r.t. the " + str(self.frame)
        return description
        
    def _new_instance(self):
        r"""
        Creates a :class:`Components` instance of the same number of indices 
        and w.r.t. the same vector frame.  

        This method must be redefined by derived classes of 
        :class:`Components`.
        
        """
        return Components(self.frame, self.nid)

    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        EXAMPLES:

        Copy on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')  
            sage: c_xy.<x,y> = m.chart('x y')
            sage: c = Components(m.default_frame(), 2)
            sage: c[:] = [[1, x+y], [x*y, 2]]
            sage: c1 = c.copy() ; c1
            2-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy))
            sage: c1[:]
            [    1 x + y]
            [  x*y     2]
            sage: c1 == c
            True
            sage: c1 is c
            False
            
        The symmetries are copied::
        
            sage: c = CompFullySym(m.default_frame(), 2)
            sage: c[0,0], c[0,1], c[1,1] = (x, 2, x*y)
            sage: c[:]
            [  x   2]
            [  2 x*y]
            sage: c1 = c.copy() ; c1
            fully symmetric 2-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy))
            sage: c1[:]
            [  x   2]
            [  2 x*y]
            sage: c1 == c
            True

        """
        result = self._new_instance()
        for ind, val in self._comp.items():
             result._comp[ind] = val.copy()
        return result


    def _del_zeros(self):
        r"""
        Deletes all the zeros in the dictionary :attr:`_comp`
        
        """
        # The zeros are first searched; they are deleted in a second stage, to
        # avoid changing the dictionary while it is read
        zeros = []
        for ind, value in self._comp.items():
            if value == 0:
                zeros.append(ind)
        for ind in zeros:
            del self._comp[ind] 

    def _check_indices(self, indices):
        r"""
        Check the validity of a list of indices and returns a tuple from it
        
        INPUT:
        
        - ``indices`` -- list of indices (possibly a single integer if
          self is a 1-index object)
          
        OUTPUT:
        
        - a tuple containing valid indices
          
        """
        if isinstance(indices, (int, Integer)):
            ind = (indices,)
        else:
            ind = tuple(indices)
        if len(ind) != self.nid:
            raise TypeError("Wrong number of indices: " + str(self.nid) + 
                            " expected, while " + str(len(ind)) + 
                            " are provided.")
        n = self.manifold.dim
        si = self.manifold.sindex
        for k in range(self.nid):
            i = ind[k]
            if i < si or i > n-1+si: 
                raise IndexError("Index out of range: " 
                                  + str(i) + " not in [" + str(si) + ","
                                  + str(n-1+si) + "]")
        return ind
    
    def _get_list(self, ind_slice, chart=None):
        r"""
        Return the list of components.
        
        INPUT:
        
        - ``ind_slice`` --  a slice object
        - ``chart`` -- (default: None) chart defining the 
          coordinates used in symbolic expressions; if none is provided the 
          domain's default chart is assumed.

        
        OUTPUT:
        
        - the full list of components if  ``ind_slice`` == ``[:]``, or a slice
          of it if ``ind_slice`` == ``[a:b]`` (1-D case), in the form
          ``T[i][j]...`` for the components `T_{ij...}` (for a 2-indices 
          object, a matrix is returned).
          
        """
        from sage.matrix.constructor import matrix
        from chart import FunctionChart
        if (ind_slice.start is not None or ind_slice.stop is not None) \
          and self.nid != 1:
            raise NotImplementedError("Function [start:stop] not " +
                          "implemented for components with " + str(self.nid) + 
                          " indices.")
        si = self.manifold.sindex
        nsi = si + self.manifold.dim
        if self.nid == 1:
            if ind_slice.start is None: 
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None: 
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("Function [start:stop:step] " +
                                              "not implemented.")
            return [self[i, chart] for i in range(start, stop)]
        if self.nid == 2:
            res = [[] for i in range(si, nsi)]
            for i in range(si, nsi):
                for j in range(si, nsi):
                    val = self[i,j, chart]
                    if isinstance(val, FunctionChart):
                        # to allow for a common ring for matrix(), the
                        # symbolic expression is used instead of the 
                        # FunctionChart:
                        val = val.express  
                    res[i-si].append(val)
            return matrix(res)
        if self.nid == 3:
            return [[[ self[i,j,k, chart] for k in range(si, nsi)]
                       for j in range(si, nsi)] for i in range(si, nsi)] 
        if self.nid == 4:
            return [[[[ self[i,j,k,l, chart] for l in range(si, nsi)]
                        for k in range(si, nsi)] for j in range(si, nsi)] 
                        for i in range(si, nsi)]
        if self.nid == 5:
            return [[[[[ self[i,j,k,l,m, chart] for m in range(si, nsi)]
                         for l in range(si, nsi)] for k in range(si, nsi)] 
                         for j in range(si, nsi)] for i in range(si, nsi)]
        raise NotImplementedError("Function [:] not implemented " + 
                                       "for components with " + str(self.nid) + 
                                       " indices.")

    def _set_list(self, ind_slice, values, chart=None):
        r"""
        Set the components from a list.
        
        INPUT:
        
        - ``ind_slice`` --  a slice object
        - ``values`` -- list of values for the components : the full list if       
          ``ind_slice`` == ``[:]``, in the form ``T[i][j]...`` for the 
          component `T_{ij...}`. In the 1-D case, ``ind_slice`` can be
          a slice of the full list, in the form  ``[a:b]``
        - ``chart`` -- (default: None) chart defining the 
          coordinates used in symbolic expressions; if none is provided the 
          domain's default chart is assumed.
          
        """
        if (ind_slice.start is not None or ind_slice.stop is not None) \
          and self.nid != 1:
            raise NotImplementedError("Function [start:stop] not " +
                          "implemented for components with " + str(self.nid) + 
                          " indices.")
        if chart is None:
            chart = self.domain.def_chart
        si = self.manifold.sindex
        nsi = si + self.manifold.dim
        if self.nid == 1:
            if ind_slice.start is None: 
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None: 
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("Function [start:stop:step] " +
                                              "not implemented.")
            for i in range(start, stop):
                self[i, chart] = values[i-start]
        elif self.nid == 2:
            for i in range(si, nsi):
                for j in range(si, nsi):
                    self[i,j, chart] = values[i-si][j-si]
        elif self.nid == 3:
            for i in range(si, nsi):
                for j in range(si, nsi):
                    for k in range(si, nsi):
                        self[i,j,k, chart] = values[i-si][j-si][k-si]
        elif self.nid == 4:
            for i in range(si, nsi):
                for j in range(si, nsi):
                    for k in range(si, nsi):
                        for l in range(si, nsi):
                            self[i,j,k,l, chart] = \
                                values[i-si][j-si][k-si][l-si]
        elif self.nid == 5:
            for i in range(si, nsi):
                for j in range(si, nsi):
                    for k in range(si, nsi):
                        for l in range(si, nsi):
                            for m in range(si, nsi):
                                self[i,j,k,l,m, chart] = \
                                    values[i-si][j-si][k-si][l-si][m-si]
        else:
            raise NotImplementedError("Function [:] not implemented for " + 
                                      "components with " + str(self.nid) + 
                                      " indices.")
    
    def __getitem__(self, args):
        r"""
        Returns the component corresponding to the given indices.

        INPUT:
        
        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components.

        OUTPUT:
        
        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the components
          `T_{ij...}` (for a 2-indices object, a matrix is returned).
    
        """
        from chart import Chart
        get_scalar_fields = isinstance(args, list)
        if get_scalar_fields:
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                indices = tuple(args[0])
            else:
                indices = tuple(args)
        else:
            if isinstance(args, (int, Integer)) or isinstance(args, slice):
                chart = self.domain.def_chart
                indices = args
            elif isinstance(args[-1], Chart):
                chart = args[-1]
                indices = args[:-1]
                if len(indices) == 1:
                    indices = indices[0]
            else:
                chart = self.domain.def_chart
                indices = args
        if isinstance(indices, slice):
            return self._get_list(indices, chart)
        else: 
            # Case where indices is a set of indices
            ind = self._check_indices(indices)
            if ind in self._comp:
                if get_scalar_fields:
                    return self._comp[ind]
                else:
                    return self._comp[ind].function_chart(chart) 
            else:  # the value is zero:
                if get_scalar_fields:
                    return self.domain.zero_scalar_field
                else:
                    return chart.zero_function

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:
        
        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) ; if [:] is provided, all the components 
          are set. 
        - ``value`` -- the value to be set or a list of values if ``args``
          == ``[:]``
    
        """
        from scalarfield import ScalarField
        from chart import Chart, FunctionChart
        if isinstance(args, list):    
            # to ensure equivalence between [i,j,...] and [[i,j,...]] or 
            # [[(i,j,...)]]
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                args = tuple(args[0])
            else:
                args = tuple(args)
        if isinstance(args, (int, Integer)) or isinstance(args, slice):
            chart = self.domain.def_chart
            indices = args
        elif isinstance(args[-1], Chart):
            chart = args[-1]
            indices = args[:-1]
            if len(indices) == 1:
                indices = indices[0]
        else:
            chart = self.domain.def_chart
            indices = args
            
        if isinstance(indices, slice):
            self._set_list(indices, value, chart)
        else: 
            # Case where indices is a set of indices
            ind = self._check_indices(indices)
            if value == 0:   #!# is_zero() instead ?
                # if the component has been set previously, it is deleted,
                # otherwise nothing is done:
                if ind in self._comp:
                    del self._comp[ind]
            elif isinstance(value, ScalarField):
                self._comp[ind] = value
            elif isinstance(value, FunctionChart):
                self._comp[ind] = value.scalar_field()
            else:
                self._comp[ind] = ScalarField(self.domain, value, chart)

    def swap_adjacent_indices(self, pos1, pos2, pos3):
        r"""
        Swap two adjacent sets of indices. 
        
        This method is essentially required to reorder the covariant and 
        contravariant indices in the computation of a tensor product. 
        
        INPUT:
        
        - ``pos1`` -- position of the first index of set 1
        - ``pos2`` -- position of the first index of set 2 = 1 + position of 
          the last index of set 1 (since the two sets are adjacent)
        - ``pos3`` -- 1 + position of the last index of set 2
        
        OUTPUT:
        
        - Components with index set 1 permuted with index set 2. 
        
        EXAMPLES:
        
        Swap of the two indices of a 2-indices set of components::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = m.chart('x y z')
            sage: e = VectorFrame(m, 'e')
            sage: a = Components(e, 1)
            sage: a[:] = (1,2,3)
            sage: b = Components(e, 1)
            sage: b[:] = (4,5,6)
            sage: c = a*b ; c[:]
            [ 4  5  6]
            [ 8 10 12]
            [12 15 18]
            sage: c1 = c.swap_adjacent_indices(0,1,2) ; c1[:]
            [ 4  8 12]
            [ 5 10 15]
            [ 6 12 18]
            sage: c1[:] == c[:].transpose()
            True
            
        Swap of two pairs of indices on a 4-indices set of components::
            
            sage: d = c*c1 ; d 
            4-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: d1 = d.swap_adjacent_indices(0,2,4)
            sage: d[0,1,1,2]
            75
            sage: d1[1,2,0,1]
            75

        """
        result = self._new_instance()
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            result._comp[new_ind] = val 
            # the above writing is more efficient than result[new_ind] = val 
            # it does not work for the derived class CompWithSym, but for the 
            # latter, the function CompWithSym.swap_adjacent_indices will be
            # called and not the present function. 
        return result
        
    def is_zero(self):
        r""" 
        Return True if all the components are zero and False otherwise.

        EXAMPLES:
        
        A just-created set of components is initialized to zero::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = m.chart('x y z')
            sage: c = Components(m.default_frame(), 1)
            sage: c.is_zero()
            True
            sage: c[:]
            [0, 0, 0]
            sage: c[0] = 1 ; c[:]
            [1, 0, 0]
            sage: c.is_zero()
            False
            sage: c[0] = 0 ; c[:]
            [0, 0, 0]
            sage: c.is_zero()
            True
           
        It is equivalent to use the operator == to compare to zero::
        
            sage: c == 0
            True
            sage: c != 0
            False
            sage: c[0] = 1
            sage: c == 0
            False
            sage: c != 0
            True

        Comparing to a nonzero number is meaningless::
    
            sage: c == 1
            Traceback (most recent call last):
            ...
            TypeError: Cannot compare a set of components to a number.

        """
        if self._comp == {}:
            return True
        else:
            for val in self._comp.values():
                if val != 0:
                    return False
            return True

    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a set of components or 0
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if isinstance(other, (int, Integer)): # other is 0
            if other == 0:
                return self.is_zero()
            else:
                raise TypeError("Cannot compare a set of components to a " + 
                                "number.")
        else: # other is another Components
            if not isinstance(other, Components):
                raise TypeError("An instance of Components is expected.")
            if other.nid != self.nid:
                return False
            return (self - other).is_zero()

    def __ne__(self, other):
        r"""
        Inequality operator. 
        
        INPUT:
        
        - ``other`` -- a set of components or 0
        
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
        return self.copy()

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the opposite of the components represented by ``self``
    
        """
        result = self._new_instance()
        for ind, val in self._comp.items():
             result._comp[ind] = - val.copy()
        return result

    def __add__(self, other):
        r"""
        Component addition. 
        
        INPUT:
        
        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``
        
        OUTPUT:
        
        - components resulting from the addition of ``self`` and ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("The second argument for the addition must be " + 
                            "an instance of Components.")
        if isinstance(other, CompWithSym):
            return other + self     # to deal properly with the symmetries
        if other.nid != self.nid:
            raise TypeError("The two sets of components do not have the " + 
                            "same number of indices.")
        if other.frame != self.frame:
            raise TypeError("The two sets of components are not defined on " +
                            "the same vector frame.")
        result = self.copy()
        for ind, val in other._comp.items():
            result[[ind]] += val
        result._del_zeros()  #!# to be checked (maybe not necessary...) 
        return result

    def __radd__(self, other):
        r"""
        Addition on the left with ``other``. 
        
        """
        return self.__add__(other)


    def __sub__(self, other):
        r"""
        Component subtraction. 
        
        INPUT:
        
        - ``other`` -- components, of the same type as ``self``
        
        OUTPUT:
        
        - components resulting from the subtraction of ``other`` from ``self``
        
        """
        if other == 0:
            return +self
        return self.__add__(-other)  #!# correct, but probably not optimal

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``. 
        
        """
        return (-self).__add__(other)


    def __mul__(self, other):
        r"""
        Component tensor product. 
        
        INPUT:
        
        - ``other`` -- components, on the same vector frame as ``self``
        
        OUTPUT: 
        
        - the tensor product of ``self`` by ``other``
        
        """
        if not isinstance(other, Components):
            raise TypeError("The second argument for the tensor product " + 
                            "be an instance of Components.")
        if other.frame != self.frame:
            raise TypeError("The two sets of components are not defined on " +
                            "the same vector frame.")
        if isinstance(other, CompWithSym):
            sym = []
            if other.sym != []:
                for s in other.sym:
                    ns = tuple(s[i]+self.nid for i in range(len(s)))
                    sym.append(ns)
            antisym = []
            if other.antisym != []:
                for s in other.antisym:
                    ns = tuple(s[i]+self.nid for i in range(len(s)))
                    antisym.append(ns)
            result = CompWithSym(self.frame, self.nid + other.nid, sym, 
                                 antisym)
        elif self.nid == 1 and other.nid == 1:
            if self is other:  # == would be dangerous here 
                # the result is symmetric:
                result = CompFullySym(self.frame, 2)
            else:
                result = Components(self.frame, 2)
        else:
            result = Components(self.frame, self.nid + other.nid)
        for ind_s, val_s in self._comp.items():
            for ind_o, val_o in other._comp.items():
                result._comp[ind_s + ind_o] = val_s * val_o
        return result
        

    def __rmul__(self, other):
        r"""
        Multiplication on the left by ``other``. 
        
        """
        if isinstance(other, Components):
            raise NotImplementedError("Left tensor product not implemented.")
        # Left multiplication by a "scalar": 
        result = self._new_instance()
        if other == 0:
            return result   # because a Components is initialized to zero
        for ind, val in self._comp.items():
            result._comp[ind] = other * val
        result._del_zeros() #!# check (probably not necessary)
        return result


    def __div__(self, other):
        r"""
        Division (by a scalar). 
        
        """
        if isinstance(other, Components):
            raise NotImplementedError("Division by an object of type " + 
                                      "Components not implemented.")
        result = self._new_instance()
        for ind, val in self._comp.items():
            result._comp[ind] = val / other
        return result

    def mtrace(self, lpos1, lpos2):
        r"""
        Multiple index contraction (experimental). 
        The contraction is performed bet:ween corresponding 
        indices from list lpos1 and lpos2; they should be of 
        the same length. 

        INPUT:
            
        - ``lpos1`` -- first list of index positions 
        - ``lpos2`` -- second list of index positions 
          
        OUTPUT:
        
        - set of components resulting from the (lpos1, lpos2) contraction

        EXAMPLES:

        Self-contraction of a 5-indices Components::

            sage: m = Manifold(5, 'M')
            sage: coords = m.chart('x y z t u')
            sage: e = VectorFrame(m, 'e')
            sage: c = Components(e, 5)
            sage: for i in m.irange():
            ....:   for j in m.irange():
            ....:       for k in m.irange():
            ....:           for l in m.irange():
            ....:               for n in m.irange():
            ....:                   c[i,j,k,l,n] = 1
            ....:                     
            sage: f = c.mtrace([0,1], [2,3])
            sage: f[:]
            [25, 25, 25, 25, 25]        
 
        """

        # sorting - first list in natural order, second according to first
        lpos1, lpos2 = (list(i) for i in zip(*sorted(zip(lpos1, lpos2))))

        res = self 
        l = range(self.nid)

        for i in range(len(lpos1)):

            if res.nid == 2:
                a = 0 
                for j in self.manifold.irange():
                    a  += res[[j,j]]

                return a

            result = Components(self.frame, res.nid - 2)

            p1 = l.index(lpos1[i])
            p2 = l.index(lpos2[i])
    
            # reorganizing the lists
            l.remove(lpos1[i])
            l.remove(lpos2[i])

            for ind, val in res._comp.items():
                if ind[p1] == ind[p2]:
                    # there is a contribution to the contraction
                    ind_res = ind[:p1] + ind[p1+1:p2] + ind[p2+1:]
                    result[[ind_res]] += val

            res = result 
 
        return result

    def self_contract(self, pos1, pos2):
        r""" 
        Index contraction.
        
        INPUT:
            
        - ``pos1`` -- position of the first index for the contraction
        - ``pos2`` -- position of the second index for the contraction
          
        OUTPUT:
        
        - set of components resulting from the (pos1, pos2) contraction
       
        EXAMPLES:
        
        Self-contraction of a 2-indices Components::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = m.chart('x y z')
            sage: e = VectorFrame(m, 'e')
            sage: c = Components(e, 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]] 
            sage: c.self_contract(0,1).expr()
            15
            sage: c[0,0] + c[1,1] + c[2,2]  # check
            15

        Self-contraction of a 3-indices Components::
        
            sage: a = Components(e, 1)
            sage: a[:] = (1,2,3)
            sage: b = a*c ; b
            3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s = b.self_contract(0,1) ; s
            1-index components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]   
            [30, 36, 42]
            sage: [sum(b[j,j,i] for j in m.irange()) for i in m.irange()]  # check
            [30, 36, 42]            
            sage: t = b.self_contract(0,2) ; t[:]
            [14, 32, 50]
            sage: [sum(b[j,i,j] for j in m.irange()) for i in m.irange()]  # check
            [14, 32, 50]

        """
        if self.nid < 2:
            raise TypeError("Contraction can be perfomed only on " + 
                                "components with at least 2 indices.")
        if pos1 < 0 or pos1 > self.nid - 1:
            raise IndexError("pos1 out of range.")
        if pos2 < 0 or pos2 > self.nid - 1:
            raise IndexError("pos2 out of range.")
        if pos1 == pos2:
            raise IndexError("The two positions must differ for the " +
                                 "contraction to be meaningful.")
        si = self.manifold.sindex
        nsi = si + self.manifold.dim
        if self.nid == 2:
            res = 0 
            for i in range(si, nsi):
                res += self[[i,i]]
            return res
        else:
            # More than 2 indices
            result = Components(self.frame, self.nid - 2)
            if pos1 > pos2:
                pos1, pos2 = (pos2, pos1)
            for ind, val in self._comp.items():
                if ind[pos1] == ind[pos2]:
                    # there is a contribution to the contraction
                    ind_res = ind[:pos1] + ind[pos1+1:pos2] + ind[pos2+1:]
                    result[[ind_res]] += val
            result._del_zeros()  #!# maybe not necessary
            return result


    def contract(self, pos1, other, pos2):
        r""" 
        Index contraction with another instance of :class:`Components`. 
        
        INPUT:
            
        - ``pos1`` -- position of the first index (in ``self``) for the 
          contraction
        - ``other`` -- the set of components to contract with
        - ``pos2`` -- position of the second index (in ``other``) for the 
          contraction
          
        OUTPUT:
        
        - set of components resulting from the (pos1, pos2) contraction
       
        EXAMPLES:

        Contraction of a 1-index Components with a 2-index one::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = m.chart('x y z')
            sage: e = VectorFrame(m, 'e')
            sage: c = Components(e, 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]] 
            sage: a = Components(e, 1)
            sage: a[:] = (1,2,3)
            sage: s = a.contract(0,c,0) ; s
            1-index components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]
            [30, 36, 42]
            sage: [sum(a[j]*c[j,i] for j in m.irange()) for i in m.irange()]  # check
            [30, 36, 42]
            sage: t = a.contract(0,c,1) ; t[:]
            [14, 32, 50]
            sage: [sum(a[j]*c[i,j] for j in m.irange()) for i in m.irange()]  # check
            [14, 32, 50]

        """
        if not isinstance(other, Components):
            raise TypeError("For the contraction, other must be an instance " +
                            "of Components.")
        if pos1 < 0 or pos1 > self.nid - 1:
            raise IndexError("pos1 out of range.")
        if pos2 < 0 or pos2 > other.nid - 1:
            raise IndexError("pos2 out of range.")
        return (self*other).self_contract(pos1, 
                                          pos2+self.nid) #!# correct but not optimal

    def non_redundant_index_generator(self):
        r"""
        Generator of indices. 
        
        In the absence of declared symmetries, all possible indices are 
        generated. 
        
        Only versions of this method for derived classes with symmetries or
        antisymmetries are not trivial. For the base class 
        :class:`Components` (i.e. with no symmetry declated), the output is 
        simply that generated by :meth:`Manifold.index_generator`. 
                
        OUTPUT:
        
        - an iterable index

        EXAMPLES:
        
        Indices on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = m.chart('x y z')
            sage: c = Components(m.default_frame(), 1)
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1,)
            (2,)
            (3,)
            sage: c = Components(m.default_frame(), 2)
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 1)
            (1, 2)
            (1, 3)
            (2, 1)
            (2, 2)
            (2, 3)
            (3, 1)
            (3, 2)
            (3, 3)
  
        """
        # In the present case, we simply rely on Manifold.index_generator:
        for ind in self.manifold.index_generator(self.nid):
            yield ind

    def symmetrize(self, pos=None):
        r"""
        Symmetrization over the given index positions
        
        INPUT:
        
        - ``pos`` -- (default: None) list of index positions involved in the 
          symmetrization (with the convention position=0 for the first index); 
          if none, the symmetrization is performed over all the indices
          
        OUTPUT:
        
        - an instance of :class:`CompWithSym` describing the symmetrized 
          components. 
          
        EXAMPLES:
        
        Symmetrization of 2-indices components on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = m.chart('x y z')  
            sage: c = Components(m.default_frame(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s = c.symmetrize() ; s
            fully symmetric 2-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: c[:], s[:]
            (
            [1 2 3]  [1 3 5]
            [4 5 6]  [3 5 7]
            [7 8 9], [5 7 9]
            )
            sage: c.symmetrize() == c.symmetrize([0,1]) # since in the present case, there are only two indices
            True
            
        Symmetrization of 3-indices components::
        
            sage: c = Components(m.default_frame(), 3)
            sage: c[:] = [[[1,2,3], [4,5,6], [7,8,9]], [[10,11,12], [13,14,15], [16,17,18]], [[19,20,21], [22,23,24], [25,26,27]]] 
            sage: s = c.symmetrize() ; s 
            fully symmetric 3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 16/3, 29/3], [16/3, 29/3, 14], [29/3, 14, 55/3]],
              [[16/3, 29/3, 14], [29/3, 14, 55/3], [14, 55/3, 68/3]],
              [[29/3, 14, 55/3], [14, 55/3, 68/3], [55/3, 68/3, 27]]])

        Partial symmetrization of 3-indices components::
        
            sage: s = c.symmetrize([0,1]) ; s   # symmetrization on the first two indices
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
              [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
              [[13, 14, 15], [19, 20, 21], [25, 26, 27]]])
            
            sage: s = c.symmetrize([1,2]) ; s   # symmetrization on the last two indices
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (1, 2)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 3, 5], [3, 5, 7], [5, 7, 9]],
              [[10, 12, 14], [12, 14, 16], [14, 16, 18]],
              [[19, 21, 23], [21, 23, 25], [23, 25, 27]]])

        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if pos is None:
            pos = range(self.nid)
        else:
            if len(pos) < 2:
                raise TypeError("At least two index positions must be given.")
            if len(pos) > self.nid:
                raise TypeError("Number of index positions larger than the " \
                                "total number of indices.")
        n_sym = len(pos) # number of indices involved in the symmetry
        if n_sym == self.nid:
            result = CompFullySym(self.frame, self.nid)
        else:
            result = CompWithSym(self.frame, self.nid, sym=pos)
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            #!# try/except to deal with the change list() --> domain() which
            #   occurred in Sage 5.10:
            try:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.domain())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    sum += self[[ind_perm]]
            except AttributeError:
                 for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.list())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    sum += self[[ind_perm]]               
            result[[ind]] = sum / sym_group.order()
        return result

            
    def antisymmetrize(self, pos=None):
        r"""
        Antisymmetrization over the given index positions
        
        INPUT:
        
        - ``pos`` -- (default: None) list of index positions involved in the 
          antisymmetrization (with the convention position=0 for the first 
          index); if none, the antisymmetrization is performed over all the 
          indices
          
        OUTPUT:
        
        - an instance of :class:`CompWithSym` describing the antisymmetrized 
          components. 
          
        EXAMPLES:
        
        Antisymmetrization of 2-indices components on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = m.chart('x y z')  
            sage: c = Components(m.default_frame(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s = c.antisymmetrize() ; s
            fully antisymmetric 2-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: c[:], s[:]
            (
            [1 2 3]  [ 0 -1 -2]
            [4 5 6]  [ 1  0 -1]
            [7 8 9], [ 2  1  0]
            )
            sage: c.antisymmetrize() == c.antisymmetrize([0,1])
            True
            
        Antisymmetrization of 3-indices components::
        
            sage: a = Components(m.default_frame(), 1)
            sage: a[:] = [x, y, z]
            sage: b = a*c ; b  # tensor product of a and c
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s = b.antisymmetrize() ; s
            fully antisymmetric 3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s == b.antisymmetrize([0,1,2])  # s results from the antisymmetrization over all indices
            True
            sage: b[:], s[:]
            ([[[x, 2*x, 3*x], [4*x, 5*x, 6*x], [7*x, 8*x, 9*x]],
              [[y, 2*y, 3*y], [4*y, 5*y, 6*y], [7*y, 8*y, 9*y]],
              [[z, 2*z, 3*z], [4*z, 5*z, 6*z], [7*z, 8*z, 9*z]]],
             [[[0, 0, 0], [0, 0, -1/3*x + 2/3*y - 1/3*z], [0, 1/3*x - 2/3*y + 1/3*z, 0]],
              [[0, 0, 1/3*x - 2/3*y + 1/3*z], [0, 0, 0], [-1/3*x + 2/3*y - 1/3*z, 0, 0]],
              [[0, -1/3*x + 2/3*y - 1/3*z, 0], [1/3*x - 2/3*y + 1/3*z, 0, 0], [0, 0, 0]]])
            sage: s[1,2,3]
            -1/3*x + 2/3*y - 1/3*z
            sage: s[1,2,3] == (b[1,2,3]-b[1,3,2]+b[2,3,1]-b[2,1,3]+b[3,1,2]-b[3,2,1])/6
            True

        Partial antisymmetrization of 3-indices components::
            
            sage: s = b.antisymmetrize([0,1]) ; s   # antisymmetrization on the first two indices
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (0, 1)
            sage: b[:], s[:]
            ([[[x, 2*x, 3*x], [4*x, 5*x, 6*x], [7*x, 8*x, 9*x]],
              [[y, 2*y, 3*y], [4*y, 5*y, 6*y], [7*y, 8*y, 9*y]],
              [[z, 2*z, 3*z], [4*z, 5*z, 6*z], [7*z, 8*z, 9*z]]],
             [[[0, 0, 0],
               [2*x - 1/2*y, 5/2*x - y, 3*x - 3/2*y],
               [7/2*x - 1/2*z, 4*x - z, 9/2*x - 3/2*z]],
              [[-2*x + 1/2*y, -5/2*x + y, -3*x + 3/2*y],
               [0, 0, 0],
               [7/2*y - 2*z, 4*y - 5/2*z, 9/2*y - 3*z]],
              [[-7/2*x + 1/2*z, -4*x + z, -9/2*x + 3/2*z],
               [-7/2*y + 2*z, -4*y + 5/2*z, -9/2*y + 3*z],
               [0, 0, 0]]])
            sage: s[1,2,3]
            3*x - 3/2*y
            sage: s[1,2,3] == (b[1,2,3]-b[2,1,3])/2
            True
            
            sage: s = b.antisymmetrize([1,2]) ; s 
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (1, 2)
            sage: b[:], s[:]
            ([[[x, 2*x, 3*x], [4*x, 5*x, 6*x], [7*x, 8*x, 9*x]],
              [[y, 2*y, 3*y], [4*y, 5*y, 6*y], [7*y, 8*y, 9*y]],
              [[z, 2*z, 3*z], [4*z, 5*z, 6*z], [7*z, 8*z, 9*z]]],
             [[[0, -x, -2*x], [x, 0, -x], [2*x, x, 0]],
              [[0, -y, -2*y], [y, 0, -y], [2*y, y, 0]],
              [[0, -z, -2*z], [z, 0, -z], [2*z, z, 0]]])
            sage: s[1,2,3]
            -x
            sage: s[1,2,3] == (b[1,2,3]-b[1,3,2])/2
            True

        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if pos is None:
            pos = range(self.nid)
        else:
            if len(pos) < 2:
                raise TypeError("At least two index positions must be given.")
            if len(pos) > self.nid:
                raise TypeError("Number of index positions larger than the " \
                                "total number of indices.")
        n_sym = len(pos) # number of indices involved in the antisymmetry
        if n_sym == self.nid:
            result = CompFullyAntiSym(self.frame, self.nid)
        else:
            result = CompWithSym(self.frame, self.nid, antisym=pos)
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0 
            #!# try/except to deal with the change list() --> domain() which
            #   occurred in Sage 5.10:
            try:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.domain())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    if perm.sign() == 1:
                        sum += self[[ind_perm]]
                    else:
                        sum -= self[[ind_perm]]
            except AttributeError:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.list())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    if perm.sign() == 1:
                        sum += self[[ind_perm]]
                    else:
                        sum -= self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def tensor_field(self, tensor_type, name=None, latex_name=None):
        r"""
        Construct a tensor field that has ``self`` for components in the
        frame on which ``self`` is defined.

        INPUT:
        
        - ``tensor_type`` : the pair (k,l) defining the tensor type
        - ``name`` -- (default: None) name given to the tensor field
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:
        
        - an instance of :class:`TensorField`, with the particular subcases:
        
            * if (k,l)=(1,0), this is an instance of the subclass 
              :class:`VectorField` 
            * if (k,l)=(0,1), this is an instance of the subclass 
              :class:`OneForm`
            * if (k,l)=(1,1), this is an instance of the subclass 
              :class:`EndomorphismField`
        
        """
        from tensorfield import TensorField
        from vectorfield import VectorField
        from diffform import OneForm
        from rank2field import EndomorphismField
        k = tensor_type[0]
        l = tensor_type[1]
        if k+l != self.nid:
            raise TypeError("The tensor rank is not equal to the number of " +
                            "indices.")
        if (k,l) == (1,0):
            result = VectorField(self.domain, name=name, latex_name=latex_name)
        elif (k,l) == (0,1):
            result = OneForm(self.domain, name=name, latex_name=latex_name)
        elif (k,l) == (1,1):
            result = EndomorphismField(self.domain, name=name, 
                                                        latex_name=latex_name)
        else:
            result = TensorField(self.domain, k, l, name=name, 
                                                        latex_name=latex_name)
        result.components[self.frame] = self
        return result

            
#******************************************************************************

class CompWithSym(Components):
    r"""
    Class for storing components with respect to a given vector frame on a 
    differentiable manifold over `\RR`, taking into account symmetries or 
    antisymmetries among the indices. 
    
    The stored quantities can be tensor components or non-tensorial quantities, 
    such as connection coefficients. 
    
    Subclasses of :class:`CompWithSym` are
    
    * :class:`CompFullySym` for storing fully symmetric components.
    * :class:`CompFullyAntiSym` for storing fully antisymmetric components.

    
    INPUT:
    
    - ``frame`` -- vector frame with respect to which the components are 
      defined
    - ``nb_indices`` -- number of indices 
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      indices: each symmetry is described by a tuple containing the positions 
      of the involved indices, with the convention position=0 for the first
      index. For instance:
        * sym=(0,1) for a symmetry between the 1st and 2nd indices 
        * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
          indices and a symmetry between the 2nd, 4th and 5th indices.
    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the indices, with the same convention as for ``sym``. 
      
    EXAMPLES:

    The tensor product (operator \*) preserves the symmetries, with the 
    appropriate shift on the positions::
        
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z')
        sage: e = VectorFrame(m, 'e')
        sage: a = Components(e, 1)
        sage: a[0], a[1], a[2] = (0, 1, 2)
        sage: b = CompFullyAntiSym(e, 2)
        sage: b[0,1], b[0,2], b[1,2] = (4, 5, 6)
        sage: b[:]
        [ 0  4  5]
        [-4  0  6]
        [-5 -6  0]
        sage: c = b*a ; c
        3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2)), with antisymmetry on the index positions (0, 1)
        sage: c[:]
        [[[0, 0, 0], [0, 4, 8], [0, 5, 10]], [[0, -4, -8], [0, 0, 0], [0, 6, 12]], [[0, -5, -10], [0, -6, -12], [0, 0, 0]]]
        sage: d = a*b ; d
        3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2)), with antisymmetry on the index positions (1, 2)
        sage: d[:]
        [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 4, 5], [-4, 0, 6], [-5, -6, 0]], [[0, 8, 10], [-8, 0, 12], [-10, -12, 0]]]

    """
    def __init__(self, frame, nb_indices, sym=None, antisym=None):
        Components.__init__(self, frame, nb_indices)
        self.sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):  
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) < 2:
                    raise IndexError("At least two index positions must be " + 
                                     "provided to define a symmetry.")
                for i in isym:
                    if i<0 or i>self.nid-1:
                        raise IndexError("Invalid index position: " + str(i) +
                                         " not in [0," + str(self.nid-1) + "]")
                self.sym.append(tuple(isym))       
        self.antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):  
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) < 2:
                    raise IndexError("At least two index positions must be " + 
                                     "provided to define an antisymmetry.")
                for i in isym:
                    if i<0 or i>self.nid-1:
                        raise IndexError("Invalid index position: " + str(i) +
                                         " not in [0," + str(self.nid-1) + "]")
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
                             "index position appears more then once.")

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = str(self.nid)
        if self.nid == 1:
            description += "-index"
        else:
            description += "-indices"
        description += " components w.r.t. the " + str(self.frame)
        for isym in self.sym:
            description += ", with symmetry on the index positions " + \
                           str(tuple(isym))
        for isym in self.antisym:
            description += ", with antisymmetry on the index positions " + \
                           str(tuple(isym))
        return description
    
    def _new_instance(self):
        r"""
        Creates a :class:`CompWithSym` instance w.r.t. the same vector frame,
        and with the same number of indices and the same symmetries
        
        """
        return CompWithSym(self.frame, self.nid, self.sym, self.antisym)

    def _ordered_indices(self, indices):
        r"""
        Given a set of indices, returns a set of indices with the indices
        at the positions of symmetries or antisymmetries being ordered, 
        as well as some antisymmetry indicator.
 
        INPUT:
        
        - ``indices`` -- list of indices (possibly a single integer if
          self is a 1-index object)
        
        OUTPUT:
        
        - a pair `(s,ind)` where ind is a tuple that differs from the original 
          list of indices by a reordering at the positions of symmetries and
          antisymmetries and
            * `s=0` if the value corresponding to ``indices`` vanishes by 
              antisymmetry (repeated indices); `ind` is then set to None
            * `s=1` if the value corresponding to ``indices`` is the same as
              that corresponding to `ind`
            * `s=-1` if the value corresponding to ``indices`` is the opposite
              of that corresponding to `ind`
            
        """
        from sage.combinat.permutation import Permutation
        ind = list(self._check_indices(indices))
        for isym in self.sym:
            indsym = []
            for pos in isym:
                indsym.append(ind[pos])
            indsym_ordered = sorted(indsym)
            for k, pos in enumerate(isym):
                ind[pos] = indsym_ordered[k]

        sign = 1
        for isym in self.antisym:
            indsym = []
            for pos in isym:
                indsym.append(ind[pos])
            # Returns zero if some index appears twice:
            if len(indsym) != len(set(indsym)):
                return (0, None)
            # From here, all the indices in indsym are distinct and we need
            # to determine whether they form an even permutation of their 
            # ordered series
            indsym_ordered = sorted(indsym)
            for k, pos in enumerate(isym):
                ind[pos] = indsym_ordered[k]
            if indsym_ordered != indsym:
                # Permutation linking indsym_ordered to indsym:
                #  (the +1 is required to fulfill the convention of Permutation) 
                perm = [indsym.index(i) +1 for i in indsym_ordered]
                #c# print "indsym_ordered, indsym: ", indsym_ordered, indsym 
                #c# print "Permutation: ", Permutation(perm), " signature = ",  \
                #c#     Permutation(perm).signature()
                sign *= Permutation(perm).signature()

        ind = tuple(ind)
        return (sign, ind)

    def __getitem__(self, args):
        r"""
        Returns the component corresponding to the given indices.

        INPUT:
        
        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components.
          
        OUTPUT:
        
        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the components
          `T_{ij...}` (for a 2-indices object, a matrix is returned).
    
        """
        from chart import Chart
        get_scalar_fields = isinstance(args, list)
        if get_scalar_fields:
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                indices = tuple(args[0])
            else:
                indices = tuple(args)
        else:
            if isinstance(args, (int, Integer)) or isinstance(args, slice):
                chart = self.domain.def_chart
                indices = args
            elif isinstance(args[-1], Chart):
                chart = args[-1]
                indices = args[:-1]
                if len(indices) == 1:
                    indices = indices[0]
            else:
                chart = self.domain.def_chart
                indices = args

        if isinstance(indices, slice):
            return self._get_list(indices, chart)
        else: 
            # Case where indices is a set of indices
            sign, ind = self._ordered_indices(indices)
            if sign == 0: # the value is zero:
                if get_scalar_fields:
                    return self.domain.zero_scalar_field
                else:
                    return chart.zero_function
            else:
                if ind in self._comp:
                    if get_scalar_fields:
                        if sign == 1:
                            return self._comp[ind]
                        else: # sign = -1
                            return -self._comp[ind]
                    else:
                        if sign == 1:
                            return self._comp[ind].function_chart(chart)
                        else: # sign = -1
                            return -self._comp[ind].function_chart(chart)
                else:  # the value is zero:
                    if get_scalar_fields:
                        return self.domain.zero_scalar_field
                    else:
                        return chart.zero_function

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:
        
         - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) ; if [:] is provided, all the components 
          are set. 
        - ``value`` -- the value to be set or a list of values if ``args``
          == ``[:]``
    
        """
        from scalarfield import ScalarField
        from chart import Chart, FunctionChart
        if isinstance(args, list):    
            # to ensure equivalence between [i,j,...] and [[i,j,...]] or 
            # [[(i,j,...)]]
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                args = tuple(args[0])
            else:
                args = tuple(args)
        if isinstance(args, (int, Integer)) or isinstance(args, slice):
            chart = self.domain.def_chart
            indices = args
        elif isinstance(args[-1], Chart):
            chart = args[-1]
            indices = args[:-1]
            if len(indices) == 1:
                indices = indices[0]
        else:
            chart = self.domain.def_chart
            indices = args

        if isinstance(indices, slice):
            self._set_list(indices, value, chart)
        else: 
            # Case where indices is a set of indices
            sign, ind = self._ordered_indices(indices)
            if sign == 0:
                if value != 0:
                    raise ValueError(
                    "By antisymmetry, the component cannot have a nonzero " + 
                    "value for the indices " + str(indices))
            else:
                if value == 0:
                    # if the component has been set previously it is deleted, 
                    # otherwise nothing is done:
                    if ind in self._comp:
                        del self._comp[ind]
                elif isinstance(value, ScalarField):
                    if sign == 1:
                        self._comp[ind] = value
                    else:   # sign = -1
                        self._comp[ind] = -value
                elif isinstance(value, FunctionChart):
                    if sign == 1:
                        self._comp[ind] = value.scalar_field()
                    else:  # sign = -1
                        self._comp[ind] = (-value).scalar_field()
                else:
                    self._comp[ind] = ScalarField(self.domain, sign*value, chart)


    def swap_adjacent_indices(self, pos1, pos2, pos3):
        r"""
        Swap two adjacent sets of indices. 
        
        This method is essentially required to reorder the covariant and 
        contravariant indices in the computation of a tensor product. 
        
        The symmetries are preserved and the corresponding indices are adjusted
        consequently. 
        
        INPUT:
        
        - ``pos1`` -- position of the first index of set 1
        - ``pos2`` -- position of the first index of set 2 = 1 + position of 
          the last index of set 1 (since the two sets are adjacent)
        - ``pos3`` -- 1 + position of the last index of set 2
        
        OUTPUT:
        
        - Components with index set 1 permuted with index set 2. 
        
        EXAMPLES:
        
        Swap of the index in position 0 with the pair of indices in position 
        (1,2) in a set of components antisymmetric with respect to the indices
        in position (1,2)::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = m.chart('x y z')
            sage: c = CompWithSym(m.default_frame(), 3, antisym=(1,2))
            sage: c[0,0,1], c[0,0,2], c[0,1,2] = (1,2,3)
            sage: c[1,0,1], c[1,0,2], c[1,1,2] = (4,5,6)
            sage: c[2,0,1], c[2,0,2], c[2,1,2] = (7,8,9)
            sage: c[:]
            [[[0, 1, 2], [-1, 0, 3], [-2, -3, 0]],
            [[0, 4, 5], [-4, 0, 6], [-5, -6, 0]],
            [[0, 7, 8], [-7, 0, 9], [-8, -9, 0]]]
            sage: c1 = c.swap_adjacent_indices(0,1,3)
            sage: c1[:]
            [[[0, 0, 0], [1, 4, 7], [2, 5, 8]],
            [[-1, -4, -7], [0, 0, 0], [3, 6, 9]],
            [[-2, -5, -8], [-3, -6, -9], [0, 0, 0]]]
            sage: c.antisym
            [(1, 2)]
            sage: c1.antisym
            [(0, 1)]
            sage: c[0,1,2]
            3
            sage: c1[1,2,0]
            3
            sage: c1[2,1,0]
            -3

        """
        result = self._new_instance()
        # The symmetries:
        lpos = range(self.nid)
        new_lpos = lpos[:pos1] + lpos[pos2:pos3] + lpos[pos1:pos2] + lpos[pos3:]
        result.sym = []
        for s in self.sym:
            new_s = [new_lpos.index(pos) for pos in s]
            result.sym.append(tuple(sorted(new_s)))
        result.antisym = []
        for s in self.antisym:
            new_s = [new_lpos.index(pos) for pos in s]
            result.antisym.append(tuple(sorted(new_s)))
        # The values:
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            result[new_ind] = val  
        return result

    def __add__(self, other):
        r"""
        Component addition. 
        
        INPUT:
        
        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``
        
        OUTPUT:
        
        - components resulting from the addition of ``self`` and ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("The second argument for the addition must be a " + 
                            "an instance of Components.")
        if other.nid != self.nid:
            raise TypeError("The two sets of components do not have the " + 
                            "same number of indices.")
        if other.frame != self.frame:
            raise TypeError("The two sets of components are not defined on " +
                            "the same vector frame.")
        if isinstance(other, CompWithSym):
            # Are the symmetries of the same type ?
            diff_sym = set(self.sym).symmetric_difference(set(other.sym))
            diff_antisym = \
                set(self.antisym).symmetric_difference(set(other.antisym))
            if diff_sym == set() and diff_antisym == set():
                # The symmetries/antisymmetries are identical:
                result = self.copy()
                for ind, val in other._comp.items():
                    result[[ind]] += val
                result._del_zeros()  #!# may be unnecessary (check !)
                return result
            else:
                # The symmetries/antisymmetries are different: only the 
                # common ones are kept
                common_sym = []
                for isym in self.sym:
                    for osym in other.sym:
                        com = tuple(set(isym).intersection(set(osym)))
                        if len(com) > 1:
                            common_sym.append(com)
                common_antisym = []
                for isym in self.antisym:
                    for osym in other.antisym:
                        com = tuple(set(isym).intersection(set(osym)))
                        if len(com) > 1:
                            common_antisym.append(com)                           
                if common_sym != [] or common_antisym != []:
                    result = CompWithSym(self.frame, self.nid, common_sym, 
                                         common_antisym)
                else:
                    # no common symmetry -> the result is a generic Components:
                    result = Components(self.frame, self.nid)
        else:
            # other has no symmetry at all:
            result = Components(self.frame, self.nid)
                     
        for ind in self.manifold.index_generator(self.nid):
            result[[ind]] = self[[ind]] + other[[ind]]
        # result._del_zeros() #!# may be unnecessary (check !)
        return result


    def __mul__(self, other):
        r"""
        Component tensor product. 
        
        INPUT:
        
        - ``other`` -- components, on the same vector frame as ``self``
        
        OUTPUT: 
        
        - the tensor product of ``self`` by ``other``
        
        """
        if not isinstance(other, Components):
            raise TypeError("The second argument for the tensor product " + 
                            "be an instance of Components.")
        if other.frame != self.frame:
            raise TypeError("The two sets of components are not defined on " +
                            "the same vector frame.")
        sym = list(self.sym)
        antisym = list(self.antisym)
        
        if isinstance(other, CompWithSym):
            if other.sym != []:
                for s in other.sym:
                    ns = tuple(s[i]+self.nid for i in range(len(s)))
                    sym.append(ns)
            if other.antisym != []:
                for s in other.antisym:
                    ns = tuple(s[i]+self.nid for i in range(len(s)))
                    antisym.append(ns)
        result = CompWithSym(self.frame, self.nid + other.nid, sym, antisym)
                             
        for ind_s, val_s in self._comp.items():
            for ind_o, val_o in other._comp.items():
                result._comp[ind_s + ind_o] = val_s * val_o
        return result


    def self_contract(self, pos1, pos2):
        r""" 
        Index contraction, , taking care of the symmetries.
        
        INPUT:
            
        - ``pos1`` -- position of the first index for the contraction
        - ``pos2`` -- position of the second index for the contraction
          
        OUTPUT:
        
        - set of components resulting from the (pos1, pos2) contraction

        EXAMPLES:

        Self-contraction of symmetric 2-indices components::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = m.chart('x y z')
            sage: e = VectorFrame(m, 'e')
            sage: a = CompFullySym(e, 2)
            sage: a[:] = [[1,2,3],[2,4,5],[3,5,6]]
            sage: a.self_contract(0,1).expr()
            11

        Self-contraction of antisymmetric 2-indices components::
        
            sage: b = CompFullyAntiSym(e,2)
            sage: b[0,1], b[0,2], b[1,2] = (3, -2, 1)
            sage: b[:]
            [ 0  3 -2]
            [-3  0  1]
            [ 2 -1  0]
            sage: b.self_contract(0,1)  # must be zero by antisymmetry
            zero scalar field on the 3-dimensional manifold 'M'

        Self-contraction of 3-indices components with one symmetry::

            sage: v = Components(e,1)
            sage: v[:] = (-2, 4, -8)
            sage: c = v*b ; c
            3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2)), with antisymmetry on the index positions (1, 2)
            sage: s = c.self_contract(0,1) ; s
            1-index components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]
            [-28, 2, 8]
            sage: [sum(v[k]*b[k,i] for k in m.irange()) for i in m.irange()] # check
            [-28, 2, 8]
            sage: s = c.self_contract(1,2) ; s
            1-index components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]   # is zero by antisymmetry 
            [0, 0, 0]
            sage: c = b*v ; c
            3-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2)), with antisymmetry on the index positions (0, 1)
            sage: s = c.self_contract(0,1) ; s
            1-index components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]  # is zero by antisymmetry
            [0, 0, 0]
            sage: s = c.self_contract(1,2) ; s
            1-index components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]
            [28, -2, -8]
            sage: [sum(b[i,k]*v[k] for k in m.irange()) for i in m.irange()]  # check 
            [28, -2, -8]

        Self-contraction of 4-indices components with two symmetries::
        
            sage: c = a*b ; c
            4-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2)), with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s = c.self_contract(0,1) ; s  # the symmetry on (0,1) is lost:
            fully antisymmetric 2-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]
            [  0  33 -22]
            [-33   0  11]
            [ 22 -11   0]
            sage: s[:] == matrix( [[sum(c[k,k,i,j].expr() for k in m.irange()) for j in m.irange()] for i in m.irange()] )  # check
            True
            sage: s = c.self_contract(2,3) ; s  # the antisymmetry on (2,3) is lost:
            fully symmetric 2-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]  # the result is zero because two antisymmetric indices have been contracted
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: s = c.self_contract(1,2) ; s  # both symmetries are lost by this contraction
            2-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
            sage: s[:]
            [ 0  0  0]
            [-2  1  0]
            [-3  3 -1]
            sage: s[:] == matrix( [[sum(c[i,k,k,j].expr() for k in m.irange()) for j in m.irange()] for i in m.irange()] ) # check
            True
        
        """ 
        if self.nid < 2:
            raise TypeError("Contraction can be perfomed only on " + 
                            "components with at least 2 indices.")
        if pos1 < 0 or pos1 > self.nid - 1:
            raise IndexError("pos1 out of range.")
        if pos2 < 0 or pos2 > self.nid - 1:
            raise IndexError("pos2 out of range.")
        if pos1 == pos2:
            raise IndexError("The two positions must differ for the " +
                                "contraction to take place.")
        si = self.manifold.sindex
        nsi = si + self.manifold.dim
        if self.nid == 2:
            res = 0 
            for i in range(si, nsi):
                res += self[[i,i]]
            return res
        else:
            # More than 2 indices
            if pos1 > pos2:
                pos1, pos2 = (pos2, pos1)
            # Determination of the remaining symmetries:
            sym_res = list(self.sym)
            for isym in self.sym:
                isym_res = list(isym)
                if pos1 in isym:
                    isym_res.remove(pos1)
                if pos2 in isym: 
                    isym_res.remove(pos2)
                if len(isym_res) < 2:       # the symmetry is lost
                    sym_res.remove(isym)
                else:
                    sym_res[sym_res.index(isym)] = tuple(isym_res)
            antisym_res = list(self.antisym)
            for isym in self.antisym:
                isym_res = list(isym)
                if pos1 in isym:
                    isym_res.remove(pos1)
                if pos2 in isym: 
                    isym_res.remove(pos2)
                if len(isym_res) < 2:       # the symmetry is lost
                    antisym_res.remove(isym)
                else:
                    antisym_res[antisym_res.index(isym)] = tuple(isym_res)
            # Shift of the index positions to take into account the
            # suppression of 2 indices:
            max_sym = 0
            for k in range(len(sym_res)):
                isym_res = []
                for pos in sym_res[k]:
                    if pos < pos1:
                        isym_res.append(pos)
                    elif pos < pos2:
                        isym_res.append(pos-1)
                    else:
                        isym_res.append(pos-2)
                max_sym = max(max_sym, len(isym_res))
                sym_res[k] = tuple(isym_res)        
            max_antisym = 0
            for k in range(len(antisym_res)):
                isym_res = []
                for pos in antisym_res[k]:
                    if pos < pos1:
                        isym_res.append(pos)
                    elif pos < pos2:
                        isym_res.append(pos-1)
                    else:
                        isym_res.append(pos-2)
                max_antisym = max(max_antisym, len(isym_res))
                antisym_res[k] = tuple(isym_res)
            # Construction of the appropriate object in view of the
            # remaining symmetries:
            nid_res = self.nid - 2
            if max_sym == 0 and max_antisym == 0:
                result = Components(self.frame, nid_res)
            elif max_sym == nid_res:
                result = CompFullySym(self.frame, nid_res)
            elif max_antisym == nid_res:
                result = CompFullyAntiSym(self.frame, nid_res)
            else:
                result = CompWithSym(self.frame, nid_res,  sym=sym_res, 
                                     antisym=antisym_res)
            # The contraction itself:
            for ind_res in self.manifold.index_generator(nid_res):
                ind = list(ind_res)
                ind.insert(pos1, 0)
                ind.insert(pos2, 0)
                res = 0
                for i in range(si, nsi):
                    ind[pos1] = i
                    ind[pos2] = i 
                    res += self[[ind]]
                result[[ind_res]] = res
            result._del_zeros() #!# maybe not necessary
            return result


    def non_redundant_index_generator(self):
        r"""
        Generator of indices, with only ordered indices in case of symmetries, 
        so that only non-redundant indices are generated. 
                
        OUTPUT:
        
        - an iterable index

        EXAMPLES:
        
        Indices on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = m.chart('x y z')
            sage: c = CompFullySym(m.default_frame(), 2)
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 1)
            (1, 2)
            (1, 3)
            (2, 2)
            (2, 3)
            (3, 3)
            sage: c = CompFullyAntiSym(m.default_frame(), 2)
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 2)
            (1, 3)
            (2, 3)
            sage: c = CompWithSym(m.default_frame(), 3, sym=(1,2)) # symmetry on the last two indices
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 1, 1)
            (1, 1, 2)
            (1, 1, 3)
            (1, 2, 2)
            (1, 2, 3)
            (1, 3, 3)
            (2, 1, 1)
            (2, 1, 2)
            (2, 1, 3)
            (2, 2, 2)
            (2, 2, 3)
            (2, 3, 3)
            (3, 1, 1)
            (3, 1, 2)
            (3, 1, 3)
            (3, 2, 2)
            (3, 2, 3)
            (3, 3, 3)
            sage: c = CompWithSym(m.default_frame(), 3, antisym=(0,1)) # antisymmetry on the first two indices
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 2, 1)
            (1, 2, 2)
            (1, 2, 3)
            (1, 3, 1)
            (1, 3, 2)
            (1, 3, 3)
            (2, 3, 1)
            (2, 3, 2)
            (2, 3, 3)
            sage: c = CompFullySym(m.default_frame(), 3)
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 1, 1)
            (1, 1, 2)
            (1, 1, 3)
            (1, 2, 2)
            (1, 2, 3)
            (1, 3, 3)
            (2, 2, 2)
            (2, 2, 3)
            (2, 3, 3)
            (3, 3, 3)
            sage: c = CompFullyAntiSym(m.default_frame(), 3) 
            sage: for ind in c.non_redundant_index_generator(): print ind
            (1, 2, 3)

        """
        si = self.manifold.sindex
        imax = self.manifold.dim - 1 + si
        ind = [si for k in range(self.nid)]
        ind_end = [si for k in range(self.nid)]
        ind_end[0] = imax+1
        while ind != ind_end:
            ordered = True
            for isym in self.sym:
                for k in range(len(isym)-1):
                    if ind[isym[k+1]] < ind[isym[k]]:
                        ordered = False
                        break                
            for isym in self.antisym:
                for k in range(len(isym)-1):
                    if ind[isym[k+1]] <= ind[isym[k]]:
                        ordered = False
                        break
            if ordered:
                yield tuple(ind)
            ret = 1
            for pos in range(self.nid-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    def symmetrize(self, pos=None):
        r"""
        Symmetrization over the given index positions
        
        INPUT:
        
        - ``pos`` -- (default: None) list of index positions involved in the 
          symmetrization (with the convention position=0 for the first index); 
          if none, the symmetrization is performed over all the indices
          
        OUTPUT:
        
        - an instance of :class:`CompWithSym` describing the symmetrized 
          components. 
          
        EXAMPLES:
        
        Symmetrization of 3-indices components on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = m.chart('x y z')  
            sage: c = Components(m.default_frame(), 3)
            sage: c[:] = [[[1,2,3], [4,5,6], [7,8,9]], [[10,11,12], [13,14,15], [16,17,18]], [[19,20,21], [22,23,24], [25,26,27]]] 
            sage: cs = c.symmetrize([0,1]) ; cs  # creation of the CompWithSym object
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1)
            sage: s = cs.symmetrize() ; s   
            fully symmetric 3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: cs[:], s[:]
            ([[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
              [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
              [[13, 14, 15], [19, 20, 21], [25, 26, 27]]],
             [[[1, 16/3, 29/3], [16/3, 29/3, 14], [29/3, 14, 55/3]],
              [[16/3, 29/3, 14], [29/3, 14, 55/3], [14, 55/3, 68/3]],
              [[29/3, 14, 55/3], [14, 55/3, 68/3], [55/3, 68/3, 27]]])
            sage: s == c.symmetrize()  # should be true
            True
            sage: s1 = cs.symmetrize([0,1]) ; s1   # should return a copy of cs
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1)
            sage: s1 == cs    # check that s1 is a copy of cs
            True
        
        Let us now start with a symmetry on the last two indices::

            sage: cs1 = c.symmetrize([1,2]) ; cs1
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (1, 2)
            sage: s2 = cs1.symmetrize() ; s2
            fully symmetric 3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s2 == c.symmetrize()
            True
            
        Partial symmetrization of 4-indices components with an antisymmetry on
        the last two indices::
        
            sage: a = Components(m.default_frame(), 2)
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b = CompFullyAntiSym(m.default_frame(), 2)
            sage: b[1,2], b[1,3], b[2,3] = (2,4,8)
            sage: c = a*b ; c
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (2, 3)
            sage: s = c.symmetrize([0,1]) ; s  # symmetrization on the first two indices
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s[1,2,3,2] == (c[1,2,3,2] + c[2,1,3,2]) / 2 # check of the symmetrization 
            True
            sage: s = c.symmetrize() ; s  # symmetrization over all the indices
            fully symmetric 4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s == 0    # the full symmetrization results in zero due to the antisymmetry on the last two indices
            True
            sage: s = c.symmetrize([2,3]) ; s
            fully symmetric 4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s == 0    # must be zero since the symmetrization has been performed on the antisymmetric indices
            True
            sage: s = c.symmetrize([0,2]) ; s
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 2)
            sage: s != 0  # s is not zero, but the antisymmetry on (2,3) is lost because the position 2 is involved in the new symmetry
            True

        Partial symmetrization of 4-indices components with an antisymmetry on
        the last three indices::

            sage: a = Components(m.default_frame(), 1)
            sage: a[:] = (1, 2, 3)
            sage: b = CompFullyAntiSym(m.default_frame(), 3)
            sage: b[1,2,3] = 4
            sage: c = a*b ; c
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (1, 2, 3)
            sage: s = c.symmetrize([0,1]) ; s
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: # Note that the antisymmetry on (1, 2, 3) has been reduced to (2, 3) only
            sage: s = c.symmetrize([1,2]) ; s
            fully symmetric 4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s == 0 # because (1,2) are involved in the antisymmetry
            True

        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if pos is None:
            pos = range(self.nid)
        else:
            if len(pos) < 2:
                raise TypeError("At least two index positions must be given.")
            if len(pos) > self.nid:
                raise TypeError("Number of index positions larger than the " \
                                "total number of indices.")
        pos = tuple(pos)
        pos_set = set(pos)
        # If the symmetry is already present, there is nothing to do:
        for isym in self.sym:
            if pos_set.issubset(set(isym)):
                return self.copy()
        #
        # Interference of the new symmetry with existing ones:
        # 
        sym_res = [pos]  # starting the list of symmetries of the result
        for isym in self.sym:
            inter = pos_set.intersection(set(isym))
            # if len(inter) == len(isym), isym is included in the new symmetry
            # and therefore has not to be included in sym_res
            if len(inter) != len(isym):
                if len(inter) >= 1: 
                    # some part of isym is lost
                    isym_set = set(isym)
                    for k in inter:
                        isym_set.remove(k)
                    if len(isym_set) > 1:
                        # some part of isym remains and must be included in sym_res:
                        isym_res = tuple(isym_set)
                        sym_res.append(isym_res)
                else:
                    # case len(inter)=0: no interference: the existing symmetry is
                    # added to the list of symmetries for the result:
                    sym_res.append(isym)
        #
        # Interference of the new symmetry with existing antisymmetries:
        # 
        antisym_res = []  # list of antisymmetries of the result
        for iasym in self.antisym:
            inter = pos_set.intersection(set(iasym))
            if len(inter) > 1: 
                # If at least two of the symmetry indices are already involved 
                # in the antisymmetry, the outcome is zero: 
                return CompFullySym(self.frame, self.nid)
                # (Note that a new instance of CompFullySym is initialized to zero)
            elif len(inter) == 1:
                # some piece of antisymmetry is lost
                k = inter.pop()  # the symmetry index position involved in the antisymmetry
                iasym_set = set(iasym)
                iasym_set.remove(k)
                if len(iasym_set) > 1:
                    iasym_res = tuple(iasym_set)
                    antisym_res.append(iasym_res)
                # len(iasym_set) = 1, the antisymmetry is fully lost, it is 
                # therefore not appended to antisym_res
            else:
                # case len(inter)=0: no interference: the antisymmetry is
                # added to the list of antisymmetries for the result:
                antisym_res.append(iasym)
        #
        # Creation of the result object
        #
        max_sym = 0
        for isym in sym_res:
            max_sym = max(max_sym, len(isym))
        if max_sym == self.nid:
            result = CompFullySym(self.frame, self.nid)
        else:
            result = CompWithSym(self.frame, self.nid, sym=sym_res, 
                                 antisym=antisym_res)
        #
        # Symmetrization
        #
        n_sym = len(pos) # number of indices involved in the symmetry
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0 
            #!# try/except to deal with the change list() --> domain() which
            #   occurred in Sage 5.10:
            try:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.domain())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    sum += self[[ind_perm]]
            except AttributeError:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.list())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    sum += self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result


    def antisymmetrize(self, pos=None):
        r"""
        Antisymmetrization over the given index positions
        
        INPUT:
        
        - ``pos`` -- (default: None) list of index positions involved in the 
          antisymmetrization (with the convention position=0 for the first index); 
          if none, the antisymmetrization is performed over all the indices
          
        OUTPUT:
        
        - an instance of :class:`CompWithSym` describing the antisymmetrized 
          components. 
          
        EXAMPLES:
        
        Antisymmetrization of 3-indices components on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = m.chart('x y z')
            sage: a = Components(m.default_frame(), 1)
            sage: a[:] = [x, y, z]
            sage: b = CompFullyAntiSym(m.default_frame(), 2)
            sage: b[1,2], b[1,3], b[2,3] = (1,2,3)
            sage: c = a*b ; c   # tensor product of a by b
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (1, 2)
            sage: s = c.antisymmetrize() ; s     
            fully antisymmetric 3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: c[:], s[:]
            ([[[0, x, 2*x], [-x, 0, 3*x], [-2*x, -3*x, 0]],
              [[0, y, 2*y], [-y, 0, 3*y], [-2*y, -3*y, 0]],
              [[0, z, 2*z], [-z, 0, 3*z], [-2*z, -3*z, 0]]],
             [[[0, 0, 0], [0, 0, x - 2/3*y + 1/3*z], [0, -x + 2/3*y - 1/3*z, 0]],
              [[0, 0, -x + 2/3*y - 1/3*z], [0, 0, 0], [x - 2/3*y + 1/3*z, 0, 0]],
              [[0, x - 2/3*y + 1/3*z, 0], [-x + 2/3*y - 1/3*z, 0, 0], [0, 0, 0]]])
            sage: s[1,2,3]
            x - 2/3*y + 1/3*z
            sage: s[1,2,3] == (c[1,2,3]-c[1,3,2]+c[2,3,1]-c[2,1,3]+c[3,1,2]-c[3,2,1])/6
            True
            
        Partial antisymmetrization: the antisymmetrization on the first two indices erases the 
        existing antisymmetry on the last two ones::
        
            sage: s = c.antisymmetrize([0,1]) ; s  
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (0, 1)
            sage: c[:], s[:]
            ([[[0, x, 2*x], [-x, 0, 3*x], [-2*x, -3*x, 0]],
              [[0, y, 2*y], [-y, 0, 3*y], [-2*y, -3*y, 0]],
              [[0, z, 2*z], [-z, 0, 3*z], [-2*z, -3*z, 0]]],
             [[[0, 0, 0], [-1/2*x, -1/2*y, 3/2*x - y], [-x, -3/2*x - 1/2*z, -z]],
              [[1/2*x, 1/2*y, -3/2*x + y], [0, 0, 0], [-y + 1/2*z, -3/2*y, -3/2*z]],
              [[x, 3/2*x + 1/2*z, z], [y - 1/2*z, 3/2*y, 3/2*z], [0, 0, 0]]])
            sage: s[1,2,3]
            3/2*x - y
            sage: s[1,2,3] == (c[1,2,3]-c[2,1,3])/2
            True
            sage: s = c.antisymmetrize([1,2]) ; s  # antisymmetrization on the last two indices
            3-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (1, 2)
            sage: s == c
            True

        Partial antisymmetrization of 4-indices components with a symmetry on 
        the first two indices::
            
            sage: a = CompFullySym(m.default_frame(), 2)
            sage: a[1,1], a[1,2], a[1,3] = (x, y, z) 
            sage: a[2,2], a[2,3], a[3,3] = (y^2, z^2, z^3)
            sage: b = Components(m.default_frame(), 2)
            sage: b[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: c = a*b ; c
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1)
            sage: s = c.antisymmetrize([2,3]) ; s
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s[3,3,1,2]
            -z^3
            sage: s[3,3,1,2] == (c[3,3,1,2]-c[3,3,2,1])/2
            True

        The full antisymmetrization results in zero because of the symmetry on the 
        first two indices::
        
            sage: s = c.antisymmetrize() ; s
            fully antisymmetric 4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s == 0
            True
            
        Similarly, the partial antisymmetrization on the first two indices results in zero::
        
            sage: s = c.antisymmetrize([0,1]) ; s
            fully antisymmetric 4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: s == 0
            True
            
        The partial antisymmetrization on the positions (0,2) destroys the symmetry on (0,1)::
        
            sage: s = c.antisymmetrize([0,2]) ; s
            4-indices components w.r.t. the coordinate frame (M, (d/dx,d/dy,d/dz)), with antisymmetry on the index positions (0, 2)
            sage: s != 0
            True
            sage: s[1,2,3,1]
            -1/2*z^2 + 7/2*y
            sage: s[2,1,3,1]  # the symmetry (0,1) is lost
            7/2*y - 2*z
            sage: s[3,2,1,1]  # the antisymmetry (0,2) holds
            1/2*z^2 - 7/2*y

        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if pos is None:
            pos = range(self.nid)
        else:
            if len(pos) < 2:
                raise TypeError("At least two index positions must be given.")
            if len(pos) > self.nid:
                raise TypeError("Number of index positions larger than the " \
                                "total number of indices.")
        pos = tuple(pos)
        pos_set = set(pos)
        # If the antisymmetry is already present, there is nothing to do:
        for iasym in self.antisym:
            if pos_set.issubset(set(iasym)):
                return self.copy()
        #
        # Interference of the new antisymmetry with existing ones
        # 
        antisym_res = [pos]  # starting the list of symmetries of the result
        for iasym in self.antisym:
            inter = pos_set.intersection(set(iasym))
            # if len(inter) == len(iasym), iasym is included in the new 
            # antisymmetry and therefore has not to be included in antisym_res
            if len(inter) != len(iasym):
                if len(inter) >= 1: 
                    # some part of iasym is lost
                    iasym_set = set(iasym)
                    for k in inter:
                        iasym_set.remove(k)
                    if len(iasym_set) > 1:
                        # some part of iasym remains and must be included in 
                        # antisym_res:
                        iasym_res = tuple(iasym_set)
                        antisym_res.append(iasym_res)
                else:
                    # case len(inter)=0: no interference: the existing 
                    # antisymmetry is added to the list of antisymmetries for 
                    # the result:
                    antisym_res.append(iasym)
        #
        # Interference of the new antisymmetry with existing symmetries:
        # 
        sym_res = []  # list of symmetries of the result
        for isym in self.sym:
            inter = pos_set.intersection(set(isym))
            if len(inter) > 1: 
                # If at least two of the antisymmetry indices are already involved 
                # in the symmetry, the outcome is zero: 
                return CompFullyAntiSym(self.frame, self.nid)
                # (Note that a new instance of CompFullyAntiSym is initialized to zero)
            elif len(inter) == 1:
                # some piece of the symmetry is lost
                k = inter.pop()  # the antisymmetry index position involved in the symmetry
                isym_set = set(isym)
                isym_set.remove(k)
                if len(isym_set) > 1:
                    isym_res = tuple(isym_set)
                    sym_res.append(isym_res)
                # len(isym_set) = 1, the symmetry is fully lost, it is 
                # therefore not appended to sym_res
            else:
                # case len(inter)=0: no interference: the symmetry is
                # added to the list of symmetries for the result:
                sym_res.append(isym)
        #
        # Creation of the result object
        #
        max_sym = 0
        for isym in antisym_res:
            max_sym = max(max_sym, len(isym))
        if max_sym == self.nid:
            result = CompFullyAntiSym(self.frame, self.nid)
        else:
            result = CompWithSym(self.frame, self.nid, sym=sym_res, 
                                 antisym=antisym_res)
        #
        # Antisymmetrization
        #
        n_sym = len(pos) # number of indices involved in the antisymmetry
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0 
            #!# try/except to deal with the change list() --> domain() which
            #   occurred in Sage 5.10:
            try:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.domain())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    if perm.sign() == 1:
                        sum += self[[ind_perm]]
                    else:
                        sum -= self[[ind_perm]]
            except AttributeError:
                for perm in sym_group.list():
                    # action of the permutation on [0,1,...,n_sym-1]:
                    perm_action = map(lambda x: x-1, perm.list())
                    ind_perm = list(ind)
                    for k in range(n_sym):
                        ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                    if perm.sign() == 1:
                        sum += self[[ind_perm]]
                    else:
                        sum -= self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def tensor_field(self, tensor_type, name=None, latex_name=None):
        r"""
        Construct a tensor field that has ``self`` for components in the
        frame on which ``self`` is defined.

        INPUT:
        
        - ``tensor_type`` : the pair (k,l) defining the tensor type
        - ``name`` -- (default: None) name given to the tensor field
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:
        
        - an instance of :class:`TensorField`
        
        """
        from tensorfield import TensorField
        k = tensor_type[0]
        l = tensor_type[1]
        if k+l != self.nid:
            raise TypeError("The tensor rank is not equal to the number of " +
                            "indices.")
        result = TensorField(self.domain, k, l, name=name, 
                             latex_name=latex_name, sym=self.sym, 
                             antisym=self.antisym)
        result.components[self.frame] = self
        return result


#******************************************************************************

class CompFullySym(CompWithSym):
    r"""
    Class for storing fully symmetric components with respect to a given vector
    frame on a differentiable manifold over `\RR`.
    
    The stored quantities can be tensor components or non-tensorial quantities, 
    such as connection coefficients. 
    
    INPUT:
    
    - ``frame`` -- vector frame with respect to which the components are 
      defined
    - ``nb_indices`` -- number of indices 
      
    """
    def __init__(self, frame, nb_indices):
        CompWithSym.__init__(self, frame, nb_indices, sym=range(nb_indices))

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return "fully symmetric " + str(self.nid) + "-indices" + \
              " components w.r.t. the " + str(self.frame)
    
    def _new_instance(self):
        r"""
        Creates a :class:`CompFullySym` instance w.r.t. the same vector frame,
        and with the same number of indices.
        
        """
        return CompFullySym(self.frame, self.nid)

    def __getitem__(self, args):
        r"""
        Returns the component corresponding to the given indices.

        INPUT:
        
        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components.

          
        OUTPUT:
        
        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the components
          `T_{ij...}` (for a 2-indices object, a matrix is returned).
    
        """
        from chart import Chart
        get_scalar_fields = isinstance(args, list)
        if get_scalar_fields:
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                indices = tuple(args[0])
            else:
                indices = tuple(args)
        else:
            if isinstance(args, (int, Integer)) or isinstance(args, slice):
                chart = self.domain.def_chart
                indices = args
            elif isinstance(args[-1], Chart):
                chart = args[-1]
                indices = args[:-1]
                if len(indices) == 1:
                    indices = indices[0]
            else:
                chart = self.domain.def_chart
                indices = args

        if isinstance(indices, slice):
            return self._get_list(indices, chart)
        else: 
            # Case where indices is a set of indices
            ind = self._ordered_indices(indices)[1]  # [0]=sign is not used
            if ind in self._comp:
                if get_scalar_fields:
                    return self._comp[ind]
                else:
                    return self._comp[ind].function_chart(chart)
            else: # the value is zero:
                if get_scalar_fields:
                    return self.domain.zero_scalar_field
                else:
                    return chart.zero_function

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:
        
        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) ; if [:] is provided, all the components 
          are set. 
        - ``value`` -- the value to be set or a list of values if ``args``
          == ``[:]``
    
        """
        from scalarfield import ScalarField
        from chart import Chart, FunctionChart
        if isinstance(args, list):    
            # to ensure equivalence between [i,j,...] and [[i,j,...]] or 
            # [[(i,j,...)]]
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                args = tuple(args[0])
            else:
                args = tuple(args)
        if isinstance(args, (int, Integer)) or isinstance(args, slice):
            chart = self.domain.def_chart
            indices = args
        elif isinstance(args[-1], Chart):
            chart = args[-1]
            indices = args[:-1]
            if len(indices) == 1:
                indices = indices[0]
        else:
            chart = self.domain.def_chart
            indices = args
            
        if isinstance(indices, slice):
            self._set_list(indices, value, chart)
        else: 
            # Case where indices is a set of indices
            ind = self._ordered_indices(indices)[1]  # [0]=sign is not used
            if value == 0:
                # if the component has been set previously it is deleted, 
                # otherwise nothing is done: 
                if ind in self._comp:
                    del self._comp[ind]
            elif isinstance(value, ScalarField):
                self._comp[ind] = value
            elif isinstance(value, FunctionChart):
                self._comp[ind] = value.scalar_field()
            else:
                self._comp[ind] = ScalarField(self.domain, value, chart)


    def __add__(self, other):
        r"""
        Component addition. 
        
        INPUT:
        
        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``
        
        OUTPUT:
        
        - components resulting from the addition of ``self`` and ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("The second argument for the addition must be a " + 
                            "an instance of Components.")
        if isinstance(other, CompFullySym):
            if other.nid != self.nid:
                raise TypeError("The two sets of components do not have the " + 
                                "same number of indices.")
            if other.frame != self.frame:
                raise TypeError("The two sets of components are not defined " +
                                "on the same vector frame.")
            result = self.copy()
            for ind, val in other._comp.items():
                result[[ind]] += val
            result._del_zeros()  #!# may be unnecessary (check !)
            return result
        else:
            return CompWithSym.__add__(self, other)

    def tensor_field(self, tensor_type, name=None, latex_name=None):
        r"""
        Construct a tensor field that has ``self`` for components in the
        frame on which ``self`` is defined.

        INPUT:
        
        - ``tensor_type`` : the pair (k,l) defining the tensor type
        - ``name`` -- (default: None) name given to the tensor field
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:
        
        - an instance of :class:`TensorField`; if (k,l)=(0,2), this is actually
          an instance of the subclass :class:`SymBilinFormField`
        
        """
        from tensorfield import TensorField
        from rank2field import SymBilinFormField
        k = tensor_type[0]
        l = tensor_type[1]
        if k+l != self.nid:
            raise TypeError("The tensor rank is not equal to the number of " +
                            "indices.")
        if (k,l) == (0,2):
            result = SymBilinFormField(self.domain, name=name, 
                                                        latex_name=latex_name)
        else:
            result = TensorField(self.domain, k, l, name=name, 
                                 latex_name=latex_name, sym=self.sym, 
                                 antisym=self.antisym)
        result.components[self.frame] = self
        return result


#******************************************************************************

class CompFullyAntiSym(CompWithSym):
    r"""
    Class for storing fully antisymmetric components with respect to a given 
    vector frame on a differentiable manifold over `\RR`.
    
    The stored quantities can be tensor components or non-tensorial quantities, 
    such as connection coefficients. 
    
    INPUT:
    
    - ``frame`` -- vector frame with respect to which the components are 
      defined
    - ``nb_indices`` -- number of indices 
      
    """
    def __init__(self, frame, nb_indices):
        CompWithSym.__init__(self, frame, nb_indices, 
                             antisym=range(nb_indices))

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return "fully antisymmetric " + str(self.nid) + "-indices" + \
               " components w.r.t. the " + str(self.frame)
    
    def _new_instance(self):
        r"""
        Creates a :class:`CompFullyAntiSym` instance w.r.t. the same vector frame,
        and with the same number of indices.
        
        """
        return CompFullyAntiSym(self.frame, self.nid)


    def __add__(self, other):
        r"""
        Component addition. 
        
        INPUT:
        
        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``
        
        OUTPUT:
        
        - components resulting from the addition of ``self`` and ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("The second argument for the addition must be a " + 
                            "an instance of Components.")
        if isinstance(other, CompFullyAntiSym):
            if other.nid != self.nid:
                raise TypeError("The two sets of components do not have the " + 
                                "same number of indices.")
            if other.frame != self.frame:
                raise TypeError("The two sets of components are not defined " +
                                "on the same vector frame.")
            result = self.copy()
            for ind, val in other._comp.items():
                result[[ind]] += val
            result._del_zeros()  #!# may be unnecessary (check !)
            return result
        else:
            return CompWithSym.__add__(self, other)


    def tensor_field(self, tensor_type, name=None, latex_name=None):
        r"""
        Construct a tensor field that has ``self`` for components in the
        frame on which ``self`` is defined.

        INPUT:
        
        - ``tensor_type`` : the pair (k,l) defining the tensor type
        - ``name`` -- (default: None) name given to the tensor field
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:
        
        - an instance of :class:`TensorField`; if k=0, this is actually
          an instance of the subclass :class:`DiffForm`
        
        """
        from tensorfield import TensorField
        from diffform import DiffForm
        k = tensor_type[0]
        l = tensor_type[1]
        if k+l != self.nid:
            raise TypeError("The tensor rank is not equal to the number of " +
                            "indices.")
        if k == 0:
            result = DiffForm(self.domain, self.nid, name=name, 
                                                        latex_name=latex_name)
        else:
            result = TensorField(self.domain, k, l, name=name, 
                                 latex_name=latex_name, sym=self.sym, 
                                 antisym=self.antisym)
        result.components[self.frame] = self
        return result

            
#******************************************************************************

class KroneckerDelta(CompFullySym):
    r"""
    Kronecker delta `\delta_{ij}`.
    
        
    INPUT:
    
    - ``frame`` -- vector frame with respect to which the components are 
      defined

    EXAMPLES:

    The Kronecker delta on a 3-dimensional manifold::
        
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z')
        sage: d = KroneckerDelta(m.default_frame()) ; d
        Kronecker delta of size 3x3
        sage: d[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: isinstance(d, CompFullySym) # The Kronecker delta is (fully) symmetric
        True
        
    One can read, but not set, the components of a Kronecker delta::
    
        sage: d[1,1]
        1
        sage: d[1,1] = 2
        Traceback (most recent call last):
        ...
        NotImplementedError: The components of a Kronecker delta cannot be changed.        

    """
    def __init__(self, frame):
        CompFullySym.__init__(self, frame, 2)
        from scalarfield import ScalarField
        chart = self.domain.def_chart
        for i in self.manifold.irange():
            self._comp[(i,i)] = ScalarField(self.domain, 1, chart)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        n = str(self.manifold.dim)
        return "Kronecker delta of size " + n + "x" + n  
    
    def __setitem__(self, args, value):
        r"""
        Should not be used (the components of a Kronecker delta are constant)
        """
        raise NotImplementedError("The components of a Kronecker delta " + 
                                  "cannot be changed.")
