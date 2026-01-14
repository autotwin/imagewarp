"""
ImageWarp CLI Application Runner

This module sets up and runs the command-line interface (CLI) for
the ImageWarp toolkit. It employs Typer along with custom components
(such as CustomCLIGroup and Color)  to provide enhanced, user-friendly
command processing, including rich help messages, error handling, and
color-coded outputs.

Subcommands Included:
    - hello_world: Provides minimum working examples for demonstration.

Usage:
    This module is designed to be the primary entry point from the
    command line:
    ```bash
    warp
    ```
"""

import typer
import subprocess
from pathlib import Path

from imagewarp.cli.customization import CustomCLIGroup, Color

import imagewarp.cli.hello_cli
import imagewarp.cli.images

app = typer.Typer(
    cls=CustomCLIGroup,
    no_args_is_help=True,
    add_completion=False,
    short_help=(f"[{Color.FRONT_MATTER}]" "ImageWarp" f"[/{Color.FRONT_MATTER}]"),
)


# Add module-level apps.
app.add_typer(
    imagewarp.cli.hello_cli.app,  # location of the module in the codebase
    name="hello_world",  # name to overwrite default typer behavior to turn "_" to "-"
    short_help="Minimum working examples for demonstration",  # short description
)
app.add_typer(
    imagewarp.cli.images.app,
    name="images",
    short_help="Generate images for userguide",
)


def main():
    """
    Entry point for the CLI application.

    This function invokes the Typer application configured with custom
    command groups and subcommands. Starts the command-line interface,
    allowing users to interact with the toolkit and utilize its
    various modules.
    """
    # No try/except here is necessary because our CustomCLIGroup.main()
    # already prints error only once.
    app(standalone_mode=False)


if __name__ == "__main__":
    main()
