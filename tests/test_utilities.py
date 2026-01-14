import inspect

from imagewarp import utilities as ut

from conftest import cmdrun_test


@cmdrun_test
def test_checkerboard(image_dir, request):
    n_squares = 5
    square_size_px = 10

    # get the test's own name
    fn_name = request.function.__name__
    filename = image_dir.joinpath(f"{fn_name}.png")

    board = ut.checkerboard(
        box_width=square_size_px,
        box_height=square_size_px,
        squares_x=n_squares,
        squares_y=n_squares,
    )

    ut.save_image(img=board, filename=filename)

    # assert False

    return None
