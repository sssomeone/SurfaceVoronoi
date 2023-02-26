# SurfaceVoronoi
SurfaceVoronoi: Efficiently Computing Voronoi Diagrams over Mesh Surfaces with Arbitrary Distance Solvers
Code of SurfaceVoronoi.

Abs: In this paper, we propose to compute Voronoi diagrams over mesh surfaces driven by an arbitrary geodesic distance solver, assuming that the input is a triangle mesh as well as a collection of sites $\mathbf{P}=\{p_i\}_{i=1}^m$ on the surface. We propose two key techniques to solve this problem. First, as the partition is determined by minimizing the $m$ distance fields, each of which rooted at a source site,we suggest keeping one or more distance triples, for each triangle, that may help determine the Voronoi bisectors when one uses a mark-and-sweep geodesic algorithm to predict the multi-source distance field.Second, rather than keep the distance itself at a mesh vertex, we use the squared distance to characterize the linear change of distance field restricted in a triangle, which is proved to induce an exact VD when the base surface reduces to a planar triangle mesh.Specially, our algorithm also supports the Euclidean distance, which can handle thin-sheet models (e.g. leaf) and runs fasterthan the traditional restricted Voronoi diagram~(RVD) algorithm. It is very extensible to deal with various variants of surface-based Voronoi diagrams including(1)~surface-based power diagram,(2)~constrained Voronoi diagram with curve-type breaklines,and (3)~curve-type generators. We conduct extensive experimental results to validate the ability to approximate the exact VD in different distance-driven scenarios.

Paper link: https://arxiv.org/abs/2212.09029 Doi: https://dl.acm.org/doi/abs/10.1145/3550454.3555453

## Dependence
CGAL
Eigen3
Boost

Code will coming soon.
