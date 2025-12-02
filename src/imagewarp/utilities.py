import numpy as np
from PIL import Image

def make_checkerboard(width:int, height:int, squares_x:int, squares_y: int,
                      color1=255, color2=0,):
    """
    Generate a checkerboard pattern.

    width, height : output image size in pixels
    squares_x     : number of squares horizontally
    squares_y     : number of squares vertically
    color1, color2: intensity for alternating squares (0–255 for grayscale)

    Returns a 2D uint8 NumPy array of shape (height, width).
    """
    board = np.zeros((height, width), dtype=np.uint8)
    sq_w = width  // squares_x
    sq_h = height // squares_y

    for y in range(squares_y):
        for x in range(squares_x):
            c = color1 if (x + y) % 2 == 0 else color2
            y0, y1 = y*sq_h, (y+1)*sq_h
            x0, x1 = x*sq_w, (x+1)*sq_w
            board[y0:y1, x0:x1] = c

    return board

