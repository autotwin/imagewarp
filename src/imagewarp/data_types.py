from typing import NamedTuple
from enum import Enum


# simple 2D and 3D vector types
class Vector2D(NamedTuple):
    """
    2D vector object
    """

    x: float
    y: float


class Vector3d(NamedTuple):
    """
    3D vector object
    """

    x: float
    y: float
    z: float


class TransformationType(Enum):
    """Type of transformation"""

    RIGID = (
        "rigid"  # euclidean distance preserved between points, rotation and translation
    )
    AFFINE = "affine"
