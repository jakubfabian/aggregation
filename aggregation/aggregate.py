"""
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from bisect import bisect_left
import math
import sys
from matplotlib import pyplot, colors
import numpy as np
from numpy import random
from scipy import linalg, stats
from . import generator, rotator
from .index import Index2D, Index3D


if sys.version_info[0] >= 3:
    xrange = range


def min_z_separation(elems,ref_elem,grid_res_sqr):
    """Displacement needed to connect elements.

    Compute the displacement required to move a set of spherical elements
    just below a reference element.

    Args:
        elems: The Nx3 array of element coordinates.
        ref_elem: The coordinates of the reference element.
        grid_res_sqr: The squared size of each element.

    Returns:
        The displacement.
    """
    elems = np.array(list(elems))
    if not len(elems):  
        return np.inf
    x_sep_sqr = ((elems[:,:2]-ref_elem[:2])**2).sum(1)
    match_possible = (x_sep_sqr < grid_res_sqr)
    z_sep = np.empty(match_possible.shape)
    z_sep[match_possible] = elems[match_possible,2] - ref_elem[2] - \
        np.sqrt(grid_res_sqr-x_sep_sqr[match_possible])
    z_sep[~match_possible] = np.inf            
    return z_sep.min()

def get_proj_area_from_alphashape(proj_grid, alpha=0.4):
    """
    Calculate the projected area from an alpha shape, which is based on the projected grid

    Args:
        proj_grid: 2d-projected grid
        alpha: alpha parameter: see https://pypi.org/project/alphashape/
    Returns:
        alpha shape: see https://pypi.org/project/alphashape/e

    """
    import alphashape
    from shapely.geometry import Point, Polygon
    
    coord = np.where(proj_grid>0) #get coordinates of ice pixel
 
    alpha_shape = alphashape.alphashape(np.column_stack((coord[0],coord[1])),alpha) #get alpha shape

    if proj_grid.shape[0]*proj_grid.shape[1]>10000:
        #many pixels: get area from alpha-shape polygon (computational cheaper than counting pixel and relatively accurate for large aggregates<5% deviation)
        area = alpha_shape.area 
    else: #few pixels: pixel-wise evaluation of area
        area=0
        for i in range(0,proj_grid.shape[0]):
            for j in range(0,proj_grid.shape[1]):
                if Point(i,j).intersects(alpha_shape):
                    area+=1

    return area

class Aggregate(object):
    """A volume-element aggregate model.

    This class represents a 3-D aggregate snowflake model made of many 
    volume elements.

    Constructor args:
        generator: The crystal generator used to make this aggregate.

    Constructor keyword args:
        ident: The numerical identifier for this particle 
            (integer, default 0).
    """

    def __init__(self, generator, ident=0):
        self._generator = generator
        self.grid_res = generator.grid_res            
        self.X = self._generator.generate().T
        self.ident = np.full(self.X.shape[0], ident, dtype=np.int32)        
        self.update_extent()
        self.monomer_number = 1
        self.id_tree = ident
        self.mono_type = generator.crystal


    def update_extent(self):
        """Updates the particle size information.

        This is usually handled internally and there is no need to call this
        function manually. If elements are added or removed from "X" by
        external code, this should be called after such modifications 
        are finished.
        """

        x = self.X[:,0]
        y = self.X[:,1]
        z = self.X[:,2]
        if len(x) != 0:
            self.extent = [[x.min(), x.max()], [y.min(), y.max()], 
                [z.min(), z.max()]]
        else:
            self.extent = [[0.,0.],[0.,0.],[0.,0.]]


    def project_on_dim(self, dim=2, direction=None):
        """Make a 2D projection of the aggregate.

        Args:
            dim: The dimension along which the projection is made 
                (0<=dim<=2, default 2)
            direction: a 2-tuple of Euler angles (alpha, beta) 
                in radians, giving the viewpoint direction (default None,
                if not None then dim is ignored)

        Returns:
            2D array with the projection along the given dimension.
            The projection grid spacing is equal to the aggregate
            element size. 

            If dim==0, the two dimensions of the returned array are the
            dimensions (y,z) of the aggregate (in that order).

            If dim==1, the two dimensions of the returned array are the
            dimensions (x,z) of the aggregate (in that order).

            If dim==2, the two dimensions of the returned array are the
            dimensions (x,y) of the aggregate (in that order).
        """
        if direction is not None:
            dim = 0
            (alpha, beta) = direction
            R = rotator.Rotator.rotation_matrix(alpha, beta, 0)
            self.X = self.X.dot(R)
            self.update_extent()

        try:
            ext = self.extent
            if dim == 0:
                xp = (self.X[:,1]-ext[1][0]) / self.grid_res
                yp = (self.X[:,2]-ext[2][0]) / self.grid_res
            elif dim == 1:
                xp = (self.X[:,0]-ext[0][0]) / self.grid_res
                yp = (self.X[:,2]-ext[2][0]) / self.grid_res
            elif dim == 2:
                xp = (self.X[:,0]-ext[0][0]) / self.grid_res
                yp = (self.X[:,1]-ext[1][0]) / self.grid_res
            else:
                raise AttributeError("Argument dim must be 0<=dim<=2.")

            x_max = int(round(xp.max()))
            y_max = int(round(yp.max()))

            proj_grid = np.zeros((x_max+1,y_max+1), dtype=np.uint8)
            proj_grid[xp.round().astype(int), yp.round().astype(int)] = 1

        finally:
            if direction is not None:
                self.X = self.X.dot(R.T)
                self.update_extent()

        return proj_grid


    def projected_area(self, dim=2, direction=None, method="default"):
        """Projected area of the aggregate.

        Uses the project_on_dim function to compute the projection.

        Args:
            dim: The dimension along which the projection is made 
                (0<=dim<=2, default 2)
            direction: a 2-tuple of Euler angles (alpha, beta) 
                in radians, giving the viewpoint direction (default None,
                if not None then dim is ignored)
            method: default: count pixels in projection
                    alphashape: approximate projection by an alpha-shape (this eliminates holes in the projection and smoothes the boundaries)

        Returns:
            The projected area along the given dimension.
        """
        
        if method=="default":
            proj_grid = self.project_on_dim(dim=dim, direction=direction)
            return proj_grid.sum() * self.grid_res**2
        elif method=="alphashape":

            proj_grid = self.project_on_dim(dim=dim)
            proj_grid_alpha_sum = get_proj_area_from_alphashape(proj_grid,alpha=0.5) 
            return proj_grid_alpha_sum * self.grid_res**2
        else:
            print("method \"" + method +"\" not implemented in projected area; use default()")
            proj_grid = self.project_on_dim(dim=dim, direction=direction)
            return proj_grid.sum() * self.grid_res**2
            


    def vertical_projected_area(self):
        # Deprecated, for backward compatibility
        return self.projected_area(dim=2)


    def projected_aspect_ratio(self, dim=2, direction=None):
        """The projected aspect ratio of the aggregate.

        Uses the project_on_dim function to compute the projection.

        Args:
            dim: The dimension along which the projection is made 
                (0<=dim<=2, default 2)
            direction: a 2-tuple of Euler angles (alpha, beta) 
                in radians, giving the viewpoint direction (default None,
                if not None then dim is ignored)

        Returns:
            The aspect ratio (defined as the ratio of the maximum extents
            of the projected dimensions) along the given dimension.
        """

        proj_grid = self.project_on_dim(dim=dim, direction=direction)

        x_proj = proj_grid.any(axis=0)
        y_proj = proj_grid.any(axis=1)
        x0 = np.arange(len(x_proj))[x_proj][0]
        x1 = np.arange(len(x_proj))[x_proj][-1]
        y0 = np.arange(len(y_proj))[y_proj][0]
        y1 = np.arange(len(y_proj))[y_proj][-1]
        return float(y1-y0+1)/float(x1-x0+1)

    def aspect_ratio(self):
        #calculate aspect ratio from principal_axes
        pa = self.principal_axes()
        pa_len = np.sqrt((pa**2).sum(axis=0)) # length of each axis
        width = np.sqrt(0.5*(pa_len[0]**2 + pa_len[1]**2)) #euclidean length of 2 longer eigenvectors
        height = pa_len[2]

        return height/width

    def principal_axes(self):
        """The principal axes of the aggregate.

        The principal axes are defined as the orthogonal vectors giving the
        directions of maximum variation. In other words, the aggregate can 
        be said to be largest in the direction of the first principal axis,
        and so on.

        Returns:
            A (3,3) array with a principal axis on each column, in descending
            order of length. The length of each axis gives the amount of 
            root-mean-square variation (i.e. the standard deviation) along 
            that axis.
        """

        cov = self.X.T.dot(self.X)/self.X.shape[0]
        # account for element size (this also regularizes the matrix)
        cov += np.diag(np.full(3,self.grid_res**2/12.))
        try:
            (l,v) = np.linalg.eigh(cov)
        except np.linalg.LinAlgError:
            # In case the eigenvalue computation failed (e.g. singular cov)
            v = np.zeros((3,3))
            l = np.zeros(3)
        return (v*np.sqrt(l))[:,::-1] # return in descending order

                  
    def pen_depth_intersection_mask(self, Xp, pen_depth, pen_depth_by_mass_fraction, verbose=False):
        if pen_depth_by_mass_fraction >= 100: return np.ones(len(Xp)).astype(np.bool)
        center = np.mean(Xp, axis=0)
        distance2center = np.linalg.norm(Xp - center, axis=1)
        #distance_limit = np.percentile(distance2center, pen_depth_by_mass_fraction) # gave unreliable/unreproducible results
        distance_limit = np.sort(distance2center)[int(len(distance2center)*pen_depth_by_mass_fraction/100)]
        mask_mf = distance2center <= distance_limit
        mask_pd = distance2center <= (distance2center.max() - pen_depth)
        mask = mask_mf | mask_pd
        if verbose:
            print(f"mf={sum(mask_mf)/len(mask)*100:.1f}% pd={sum(mask_pd)/len(mask)*100:.1f}% mask={sum(mask)/len(mask)*100:.1f}%")
        return mask


    def add_particle(self, particle=None, ident=None, required=False, 
        pen_depth=0.0, pen_depth_by_mass_fraction=100., add_N_monomers=None, add_id_branch=None, add_mono_type=None):

        """Merge another particle into this one.

        The other particle is added at a random location in the (x,y) plane
        and at the bottom of this particle in the z direction.

        Args:
            particle: The (N,3) array with the coordinates of the volume
                elements from the other particle.
            identifier: The (N,) array with the numerical identifiers of the
                volume elements from the other particle.
            required: Due to randomization of the merging location, a
                suitable merge point may not be found. If required==True,
                this function will keep trying until a merging point is
                found. If required==False, it will try once and then give up
                if a merging point was not found.
            pen_depth: The penetration depth, i.e. the distance that the
                other particle is allowed to penetrate inside this particle.
            pen_depth_by_mass_fraction: furthest away fraction ]1-100] of particle mass
                that is ignored when considering intersection.
                E.g. for plates, a value of .5 ignores point that are larger than radius == sqrt(.5)

        Returns:
            True if the merge was successful, False otherwise.
        """

        # measurements of the other particle
        if particle is None:
            particle = self._generator.generate().T
            add_N_monomers = 1
            add_id_branch = -9999
        x, y, z = particle.T
        extent = [(_.min(), _.max()) for _ in (x,y,z)]
        grid_res = self.grid_res

        assert pen_depth_by_mass_fraction > 1, f"{pen_depth_by_mass_fraction=} has to be between ]1,100]"
        assert pen_depth_by_mass_fraction <= 100, f"{pen_depth_by_mass_fraction=} has to be between ]1,100]"

        imask = self.pen_depth_intersection_mask(particle, pen_depth, pen_depth_by_mass_fraction)
        x, y, z = particle[imask].T
        extent = [(_.min(), _.max()) for _ in particle[imask].T]

        # limits for random positioning of the other particle
        x0 = (self.extent[0][0]-extent[0][1])
        x1 = (self.extent[0][1]-extent[0][0])
        y0 = (self.extent[1][0]-extent[1][1])
        y1 = (self.extent[1][1]-extent[1][0])        

        site_found = False
        while not site_found:
            # randomize location in x,y plane
            x_shift = x0+np.random.rand()*(x1-x0)
            y_shift = y0+np.random.rand()*(y1-y0)
            xs = x+x_shift
            ys = y+y_shift
                    
            # the overlap between this aggregate and the other particle in
            # the shifted position
            overlapping_range = \
                [max(xs.min(),self.extent[0][0])-grid_res, 
                 min(xs.max(),self.extent[0][1])+grid_res, 
                 max(ys.min(),self.extent[1][0])-grid_res, 
                 min(ys.max(),self.extent[1][1])+grid_res]
                 
            if (overlapping_range[0] >= overlapping_range[1]) or \
                (overlapping_range[2] >= overlapping_range[3]): 
                
                # no overlap, so impossible to connect -> stop
                if required:
                   continue
                else:
                   break   
        
            # elements from this particle that are candidates for connection
            imask = self.pen_depth_intersection_mask(self.X, pen_depth, pen_depth_by_mass_fraction)
            X_filter = \
                (self.X[:,0] >= overlapping_range[0]) & \
                (self.X[:,0] < overlapping_range[1]) & \
                (self.X[:,1] >= overlapping_range[2]) & \
                (self.X[:,1] < overlapping_range[3])
            overlapping_X = self.X[X_filter & imask,:]
            if not len(overlapping_X):
                if required:
                   continue
                else:
                   break
    
            # candidates from the other particle
            X_filter = \
                (xs >= overlapping_range[0]) & \
                (xs < overlapping_range[1]) & \
                (ys >= overlapping_range[2]) & \
                (ys < overlapping_range[3])
            overlapping_Xp = np.vstack((
                xs[X_filter],ys[X_filter],z[X_filter])).T
            
            # find displacement in z direction

            def find_min_distance(pc1, pc2, epsilon):
                epsilon_sqr = epsilon**2
                elem_index2d = Index2D(elem_size=epsilon)
                elem_index2d.insert(pc1[:,:2],pc1)
                min_z_sep = np.inf
                for elem in pc2:
                    # find elements in this aggregate that are near the
                    # currently tested element in the x,y plane
                    candidates = np.array(list(elem_index2d.items_near(elem[:2], epsilon)))
                    min_z_sep = min(min_z_sep, min_z_separation(candidates, elem, epsilon_sqr))
                return min_z_sep

            import os
            MAXHEIGHTCOLLISION = float(os.environ.get('MAXHEIGHTCOLLISION', 1./np.sqrt(2)))
            if MAXHEIGHTCOLLISION > 0:
                def max_height_per_bin(points, bin_size, inverse=False):
                    import numpy as np
                    """
                    Use np.digitize to bin points into a 2D grid and get the max-z point in each bin.

                    Parameters:
                        points: (N, 3) array of x, y, z
                        bin_size: scalar
                        x_range: (min_x, max_x) (optional, computed from data if None)
                        y_range: (min_y, max_y)

                    Returns:
                        (M, 3) array of highest points per bin
                    """
                    x = points[:, 0]
                    y = points[:, 1]
                    z = points[:, 2]
                    if len(points)<1: return points

                    # Set bin edges
                    x_range = (x.min(), x.max())
                    y_range = (y.min(), y.max())
                    x_bins = np.arange(x_range[0], x_range[1] + bin_size, bin_size)
                    y_bins = np.arange(y_range[0], y_range[1] + bin_size, bin_size)

                    # Bin index starts at 1 (np.digitize behavior)
                    x_idx = np.digitize(x, x_bins) - 1
                    y_idx = np.digitize(y, y_bins) - 1

                    bin_index = x_idx + y_idx * len(x_bins)

                    # Find max z per bin
                    if inverse:
                        sort_idx = np.lexsort((-z, bin_index))  # Sort by bin then z
                    else:
                        sort_idx = np.lexsort((z, bin_index))  # Sort by bin then z

                    sorted_bin_idx = bin_index[sort_idx]
                    unique, first_indices = np.unique(sorted_bin_idx[::-1], return_index=True)

                    # Get max z point per bin
                    max_indices = sort_idx[::-1][first_indices]
                    max_points = points[max_indices]

                    return max_points
                pc1 = max_height_per_bin(overlapping_X , grid_res*MAXHEIGHTCOLLISION, inverse=True)
                pc2 = max_height_per_bin(overlapping_Xp, grid_res*MAXHEIGHTCOLLISION, inverse=False)
            else:
                pc1, pc2 = overlapping_X, overlapping_Xp

            min_z_sep = find_min_distance(pc1, pc2, grid_res)

                    
            site_found = not np.isinf(min_z_sep)            
            if not required:
                break            
   
        if site_found:
            # move the candidate to the right location in the z direction
            if pen_depth_by_mass_fraction < 100:
                p_shift = particle + np.array([x_shift, y_shift, min_z_sep])
            else:
                p_shift = particle + np.array([x_shift, y_shift, min_z_sep + pen_depth])

            if ident is None:
                ident = np.zeros(p_shift.shape[0], dtype=np.int32)
            self.add_elements(p_shift, ident=ident)
            self.monomer_number += add_N_monomers
            self.id_tree = [self.id_tree, add_id_branch]
            self.mono_type = [self.mono_type, add_mono_type]
            
        return site_found        
         
   
    def align(self):
        """Align the aggregate along the principal axes.

        The longest principal axis becomes oriented along the x-axis, the 
        second longest along the y-axis, and the shortest along the z-axis.
        """

        # get and normalize principal axes
        PA = self.principal_axes()
        PA /= np.sqrt((PA**2).sum(0))
        
        # project to principal axes
        self.X = np.dot(self.X,PA)
        self.update_extent()
         
         
    def rotate(self,rotator):
        """Rotate the aggregate.

        Args:
            rotator: The rotator to be used for the rotation. See the
                rotator module.
        """
        self.X = self.X-self.X.mean(0)
        self.X = rotator.rotate(self.X.T).T
        self.update_extent()
         
    
    def visualize(self, bgcolor=(1,1,1), fgcolor=(.8,.8,.8)):
        """Visualize the aggregate using Mayavi.

        Args:
            bgcolor: Background color for the Mayavi scene.
            fgcolor: Foreground color for the Mayavi scene.
        """

        color_list = [colors.colorConverter.to_rgb(c) for c in [
            "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
            "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
            "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
        ]]

        # local import as this can take a while
        from mayavi import mlab
        
        mlab.figure(bgcolor=bgcolor, fgcolor=fgcolor)
        i = 0
        for ident in xrange(self.ident.min(), self.ident.max()+1):
            X = self.X[self.ident==ident,:]
            if X.shape[0] > 0:
                mlab.points3d(X[:,0], X[:,1], X[:,2], 
                    color=color_list[i%len(color_list)],
                    mode="cube", scale_factor=self._generator.grid_res)
                i += 1
      
      
    def grid(self, res=None):
        """Arrange elements on a regular grid.

        The gridded coordinates are the element coordinates divided by the
        res parameter and then rounded to the nearest integer. This routine
        both conserves the number of elements and gives a unique grid 
        location for each aggregate element. If more than one elements would
        end up in the same grid location, all but one are relocated into 
        nearby empty spots on the grid.

        Args:
            res: The resolution of the grid. Should be usually left at the
            default, which is the aggregate element spacing.

        Returns:
            An integer array with the gridded element coordinates as
            multiples of res.
        """


        if res==None:
            res = self.grid_res

        # This does most of the work!
        Xc = (self.X/res).round().astype(int)

        # The rest is to identify elements that would end up in the same
        # location and move them around
        N = Xc.shape[0]

        sort_ind = np.lexsort((Xc[:,2],Xc[:,1],Xc[:,0]))
        Xc = Xc[sort_ind,:]
        overlap = abs(np.diff(Xc,axis=0)).sum(1) == 0
        overlap = np.hstack((overlap, False))

        Xc_overlap = Xc[overlap,:].copy()
        Xc = Xc[~overlap,:]
        np.random.shuffle(Xc_overlap)

        for i in xrange(Xc_overlap.shape[0]):
            Xm = Xc_overlap[i,:]
            for dX in neighbors_by_distance():
                X = Xm+dX
                if not row_is_in_sorted_array(X, Xc):
                    Xc = insert_missing_row_in_sorted_array(X, Xc)
                    break

        return Xc


    def add_elements(self, added_elements, ident=0, update=True):
        """Add elements to this aggregate.

        Args:
            added_elements: A (N,3) array with the coordinates of the added
                elements.
            ident: A (N,) array with the numerical identifiers of the added
                elements.
            update: If True, the coordinates are recentered after the update.
                This is should usually be left at True, but if you call
                add_elements multiple times without calls to other Aggregate
                member functions, it will save computational effort to set
                update=False and then call update_coordinates manually
                after you're done.
        """

        self.X = np.vstack((self.X, added_elements))
        self.ident = np.hstack((self.ident, 
            np.full(added_elements.shape[0], ident, dtype=np.int32)))
        if update:
            self.update_coordinates()


    def remove_elements(self, removed_elements, tolerance=0.001, update=True):
        """Remove elements found at the given coordinates.
        
        Args:
            removed_elements: The coordinates of the elements to remove.
            tolerance: The distance from each coordinate in removed_elements,
                in multiples of grid_res, in which the elements should be 
                removed.
            update: See the update keyword argument in add_elements.
        """

        keep = np.ones(self.X.shape[0], dtype=bool)
        for re in removed_elements:
            dist_sqr = ((self.X-re)**2).sum(1)
            min_dist = dist_sqr.argmin()
            keep[dist_sqr < (self.grid_res**2 * tolerance)] = False
        self.X = self.X[keep,:]
        self.ident = self.ident[keep]
        if update:
            self.update_coordinates()


    def update_coordinates(self):
        """Recenter the aggregate and update the particle extent.
        """
        self.X -= self.X.mean(0)
        self.update_extent()


def spheres_overlap(X0, X1, r_sqr):
    return (X1[0]-X0[0])**2 + (X1[1]-X0[1])**2 + \
        (X1[2]-X0[2])**2 < r_sqr


def compare_row(x,y):
    if x[0] > y[0]:
        return 1
    elif x[0] < y[0]:
        return -1
    elif len(x)>1:
        return compare_row(x[1:],y[1:])
    else:
        return 0


def row_is_in_sorted_array(r, x):
    i0 = 0
    i1 = x.shape[0]
    while (i1-i0)>1:
        i = (i0+i1)//2
        c = compare_row(r,x[i,:])
        if c == 0:
            return True
        elif c == 1:
            i0=i
        else:
            i1=i
    return not bool(compare_row(r,x[i0,:]))


def insert_missing_row_in_sorted_array(r, x):
    i0 = 0
    i1 = x.shape[0]
    while (i1-i0)>1:
        i = (i0+i1)//2
        c = compare_row(r,x[i,:])
        if c == 1:
            i0=i
        else:
            i1=i
    if i1 > 1:
        insert_ind = i1
    else:
        insert_ind = 1 if compare_row(r,x[i0,:])==1 else 0
    
    return np.vstack((x[:insert_ind,:], r.reshape((1,3)), x[insert_ind:,:]))


def outer_layer_of_cube(cube_rad):
    """Outer layer of cube.

    Generates the coordinates on the outer layer of a cube
    centered at (0,0,0) with side of (2*cube_rad+1).
    """
    irange = list(xrange(-cube_rad, cube_rad+1))
    for dx in irange:
        for dy in irange:
            if (abs(dx)==cube_rad) or (abs(dy)==cube_rad):
                dz_iter = irange
            else:
                dz_iter = (-cube_rad, cube_rad)
            for dz in dz_iter:
                yield (dx, dy, dz)


def neighbors_by_distance():
    cube_rad = 1
    while True:
        for p in outer_layer_of_cube(cube_rad):
            yield p
        cube_rad += 1


class RimedAggregate(Aggregate):
    """A volume-element rimed aggregate model.
    
    This class adds the add_rime_particles member function to the Aggregate
    base class. See the documentation for Aggregate for more information.
    """

    RIME_IDENT = -1

    def add_rime_particles(self, N=1, pen_depth=120e-6, compact_dist=0.):
        """Add rime particles to the aggregate.
        
        Args:
            N: Number of rime particles to add.
            pen_depth: The penetration depth, i.e. the distance that the
                rime particle is allowed to penetrate inside this particle.
        """

        grid_res = self.grid_res
        grid_res_sqr = grid_res**2

        # limits for random positioning of rime particle
        x0 = (self.extent[0][0])
        x1 = (self.extent[0][1])
        y0 = (self.extent[1][0])
        y1 = (self.extent[1][1])

        use_indexing = (N > 1)

        if use_indexing:
            elem_index = Index2D(elem_size=grid_res)            
            elem_index.insert(self.X[:,:2],self.X)
            def find_overlapping(x,y,dist_mul=1):
                p_near = np.array(list(elem_index.items_near((x,y), grid_res*dist_mul)))
                if not p_near.shape[0]:
                    return p_near
                p_filter = ((p_near[:,:2]-[x,y])**2).sum(1) < grid_res_sqr*dist_mul**2
                return p_near[p_filter,:]
        else:
            def find_overlapping(x,y,dist_mul=1):
                X_filter = ((self.X[:,:2] - 
                    np.array([x,y]))**2).sum(1) < grid_res_sqr*dist_mul**2
                return self.X[X_filter,:]

        if compact_dist > 0:
            if use_indexing:
                elem_index_3d = Index3D(elem_size=grid_res)            
                elem_index_3d.insert(self.X)
                def find_overlapping_3d(x,y,z,dist_mul=1):
                    p_near = np.array(list(elem_index_3d.items_near((x,y,z), grid_res*dist_mul)))
                    if not p_near.shape[0]:
                        return p_near
                    p_filter = ((p_near-[x,y,z])**2).sum(1) < grid_res_sqr*dist_mul**2
                    return p_near[p_filter,:]
            else:
                def find_overlapping_3d(x,y,z,dist_mul=1):
                    X_filter = ((self.X - 
                        np.array([x,y,z]))**2).sum(1) < grid_res_sqr*dist_mul**2
                    return self.X[X_filter,:]

        added_particles = np.empty((N, 3))

        for particle_num in xrange(N):
            site_found = False
            while not site_found:
                xs = x0+np.random.rand()*(x1-x0)
                ys = y0+np.random.rand()*(y1-y0)

                overlapping_range = [xs-grid_res, xs+grid_res,
                     ys-grid_res, ys+grid_res]

                overlapping_X = find_overlapping(xs, ys)
                if not overlapping_X.shape[0]:
                    continue                  
                
                X_order = overlapping_X[:,2].argsort()
                overlapping_X = overlapping_X[X_order,:]            
                last_ind = bisect_left(overlapping_X[:,2], 
                    overlapping_X[0,2]+pen_depth)
                last_search_ind = bisect_left(overlapping_X[:,2], 
                    overlapping_X[0,2]+pen_depth+grid_res)
                overlapping_X = overlapping_X[:last_search_ind+1,:]
                overlapping_z = overlapping_X[:,2]

                for i in xrange(last_ind-1, -1, -1):
                    d_sqr = (overlapping_X[i,0]-xs)**2 + (overlapping_X[i,1]-ys)**2
                    dz = math.sqrt(grid_res_sqr - d_sqr)
                    z_upper = overlapping_X[i,2] + dz
                    z_lower = overlapping_X[i,2] - dz                    
                    
                    for zc in [z_upper, z_lower]:
                        overlap = False
                        if (i==0) and (zc==z_lower):
                            break # automatically attach at the last site

                        j0 = bisect_left(overlapping_z, zc-grid_res)
                        j1 = bisect_left(overlapping_z, zc+grid_res)
                        
                        #search through possible overlapping spheres
                        for j in xrange(j0, j1):
                            if j == i:
                                continue
                            elif spheres_overlap(overlapping_X[j,:], (xs,ys,zc), 
                                grid_res_sqr):
                                # there is an overlapping sphere -> site unsuitable
                                overlap = True
                                break 

                        if not overlap:
                            break

                    if not overlap:
                        # this means we found a suitable site, so add the particle

                        # run the compacting first
                        if compact_dist > 0:
                            # locate nearby particles to use for the compacting
                            X_near = find_overlapping_3d(xs, ys, zc, dist_mul=2)
                            if X_near.size > 0:
                                X = np.array([xs, ys, zc])
                                r_sqr = ((X_near-X)**2).sum(axis=1)
                                X_near = X_near[r_sqr<(2*self.grid_res)**2,:]
                                (xs, ys, zc) = self.compact_rime(X, X_near, 
                                    max_dist=compact_dist)

                        added_particles[particle_num,:] = [xs, ys, zc]
                        site_found = True
                        self.extent[0][0] = min(self.extent[0][0], xs)
                        self.extent[0][1] = max(self.extent[0][1], xs)
                        self.extent[1][0] = min(self.extent[1][0], ys)
                        self.extent[1][1] = max(self.extent[1][1], ys)
                        self.extent[2][0] = min(self.extent[2][0], zc)
                        self.extent[2][1] = max(self.extent[2][1], zc)
                        if use_indexing:
                            elem_index.insert([[xs, ys]], [[xs, ys, zc]])
                        break

        self.add_elements(added_particles, ident=self.RIME_IDENT)


    def compact_rime(self, X, X_near, max_dist=0., min_move=0.01, dr=0.1,
        max_iters=100):

        if max_dist <= 0.:
            return X

        X_old = X.copy()
        max_dist_sqr = (max_dist*self.grid_res)**2
        min_move_sqr = (min_move*self.grid_res)**2
        for it in xrange(max_iters):
            dX = X_near-X
            r_sqr = (dX**2).sum(axis=1)
            r_sqr_norm = r_sqr / self.grid_res**2
            nearest_ind = r_sqr.argsort()
            F = np.zeros(3)
            for i in nearest_ind:
                Fi = dX[i,:]/(np.sqrt(r_sqr[i])*r_sqr_norm[i])
                if r_sqr_norm[i] > 1:
                    F += Fi
                elif r_sqr_norm[i] < 0.01: #avoid singularity
                    pass
                else:
                    F -= Fi

            F *= dr
            # limit abs(F) to at most dr
            F_abs_sqr = (F**2).sum()
            if F_abs_sqr > dr**2:
                F *= dr/np.sqrt(F_abs_sqr)
            F *= self.grid_res

            X_last = X.copy()
            X += F
            dist_sqr = ((X-X_old)**2).sum()
            if dist_sqr/self.grid_res**2 > max_dist_sqr:
                # limit distance to at most max_dist
                X = X_old + (X-X_old)*(max_dist*self.grid_res/np.sqrt(dist_sqr))
                break
            if ((X-X_last)**2).sum() < min_move_sqr:
                break

        return X

         
class PseudoAggregate(Aggregate):
    """A "pseudo-aggregate" model.
        
    This is similar to the Aggregate class but instead distributes the
    crystals around by sampling the positions from a 3D normal distribution.
    See the Aggregate class for more details.
    """
    
    def __init__(self, generator, sig=1.0, ident=0):
        self.sig = sig
        self.generator = generator

        self.X = self.generator.generate().T
        x = self.X[:,0]+stats.norm.rvs(scale=sig)
        y = self.X[:,1]+stats.norm.rvs(scale=sig)
        z = self.X[:,2]+stats.norm.rvs(scale=sig)
        self.extent = [[x.min(), x.max()], [y.min(), y.max()], [z.min(), z.max()]]
        self.ident = np.full(self.X.shape[0], ident, dtype=np.int32)
        self.monomer_number = 1
        self.id_tree = ident
         
                  
    def add_particle(self, particle=None, required=False, add_N_monomers=None, add_id_branch=None):
        if particle == None:
            particle = self.generator.generate().T
            add_N_monomers = 1
        x = particle[:,0]+stats.norm.rvs(scale=self.sig)
        y = particle[:,1]+stats.norm.rvs(scale=self.sig)
        z = particle[:,2]+stats.norm.rvs(scale=self.sig)
        self.X = numpy.vstack((self.X, numpy.vstack((x,y,z)).T)) 
        x = self.X[:,0]
        y = self.X[:,1]
        z = self.X[:,2]
        self.extent = [[x.min(), x.max()], [y.min(), y.max()], [z.min(), z.max()]]
        self.monomer_number += add_N_monomers
        self.id_tree = [self.id_tree, add_id_branch]
