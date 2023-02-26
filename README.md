# SurfaceVoronoi
SurfaceVoronoi: Efficiently Computing Voronoi Diagrams over Mesh Surfaces with Arbitrary Distance Solvers
Code of SurfaceVoronoi.

Abs: Feature lines are important geometric cues in characterizing the structure of a CAD model. Despite great progress in both explicit reconstruction and implicit reconstruction, it remains a challenging task to reconstruct a polygonal surface equipped with feature lines, especially when the input point cloud is noisy and lacks faithful normal vectors. In this paper, we develop a multistage algorithm, named RFEPS, to address this challenge. The key steps include (1)denoising the point cloud based on the assumption of local planarity, (2)identifying the feature-line zone by optimization of discrete optimal transport, (3)augmenting the point set so that sufficiently many additional points are generated on potential geometry edges, and (4) generating a polygonal surface that interpolates the augmented point set based on restricted power diagram. We demonstrate through extensive experiments that RFEPS, benefiting from the edge-point augmentation and the feature-preserving explicit reconstruction, outperforms state-of-the-art methods in terms of the reconstruction quality, especially in terms of the ability to reconstruct missing feature lines.

Paper link: https://arxiv.org/abs/2212.09029 Doi: https://dl.acm.org/doi/abs/10.1145/3550454.3555453

Dependence
CGAL
Eigen3
Boost

Code will coming soon.
