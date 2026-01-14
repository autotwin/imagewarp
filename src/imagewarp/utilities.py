# default python modules

# 3rd party modules
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# # local modules
# from imagewarp.cli.customization import CustomCLIGroup, CustomCLICommand


def checkerboard(
    box_width: int,
    box_height: int,
    squares_x: int,
    squares_y: int,
    color1: int = 255,
    color2: int = 0,
) -> np.ndarray:
    """
    Generate a checkerboard.

    box_width, box_height   : size of each square (in pixels)
    squares_x, y            : number of squares horizontally and vertically
    color1, color2          : intensities for alternating squares (0–255)

    Returns a 2D uint8 array of shape
      (box_h * squares_y, box_w * squares_x).
    """
    # compute total image size
    total_w = box_width * squares_x
    total_h = box_height * squares_y

    board = np.zeros((total_h, total_w), dtype=np.uint8)

    for y in range(squares_y):
        for x in range(squares_x):
            c = color1 if (x + y) % 2 == 0 else color2
            y0, y1 = y * box_height, (y + 1) * box_height
            x0, x1 = x * box_width, (x + 1) * box_width
            board[y0:y1, x0:x1] = c

    return board


def save_image(
    img: np.ndarray,
    filename: str,
):
    plt.imshow(
        img,
        cmap="gray",
    )
    plt.axis("off")
    plt.savefig(filename, bbox_inches="tight", pad_inches=0)
    plt.close()
