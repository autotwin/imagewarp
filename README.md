# imagewarp - Homography tools for 2D images and 3D volumes

<!-- image is from here: https://www.researchgate.net/figure/Examples-of-all-the-affine-transformations-for-a-simple-geometric-shape_fig1_251509419-->

## Types of homographies (transformations)

https://www.cs.princeton.edu/courses/archive/fall99/cs426/lectures/transform/sld044.htm

https://perez-aids.medium.com/introduction-to-image-processing-part-7-homography-4af34ce9f93d

- Rigid (Euclidean): Preserves Euclidean distance between points
    - Rotations
    - Translations
    - 
- Similarity: Preserves shape of objects
    - Rigid transformations
    - Uniform scaling
- Affine: Preserves lines (alignment) and parallelism of lines
    - Similarity transformations
    - Non-uniform scaling
    - Shearing
    - Mirroring
- Projective: Preserves lines (**9 degrees of freedom**)
    - Affine transformations
    - Prespective (lines converging at a vanishing point)
    - **9 degrees of freedom**: 
