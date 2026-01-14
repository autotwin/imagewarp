"""
This module contains the CLI commands for the hello world example.
It provides simple commands to greet users and demonstrate the CLI functionality.
"""

import sys

# import click  # unused import
import typer
from imagewarp.cli.customization import CustomCLIGroup, CustomCLICommand

app = typer.Typer(
    cls=CustomCLIGroup,
    no_args_is_help=True,
    help="Hello World module – contains greeting commands.",
)

# Only force help at the module level when no subcommand is provided.
if len(sys.argv) == 2:
    sys.argv.append("--help")


@app.command(
    cls=CustomCLICommand,
    name="greet_positional",
    short_help="Greet a person with a simple hello message using positional arguments.",
)
def greet_positional(
    name: str = typer.Argument(
        ...,  # default value is empty, making it required
        help="Name of the person to greet.",
    ),
    add_on: str = typer.Argument(
        "What a lovely day!",  # default value
        help="Something else to say",
    ),
):
    """
    Greet a person with a simple hello message.

    Takes positional arguments (typer.Argument types)

    Examples for CLI usage:
        >>> warp hello_world greet_positional Andrew
        Hello, Andrew
        What a lovely day!
        >>> warp hello_world greet_positional Andrew "Clever non sequitur"
        Hello, Andrew
        Clever non sequitur

    """
    typer.echo(f"Hello, {name}!\n{add_on}")


@app.command(
    cls=CustomCLICommand,
    name="greet_keywords",
    short_help="Greet a person with a simple hello message using keyword arguments.",
)
def greet_keywords(
    name: str = typer.Option(
        ...,  # default value is empty, making it required
        help="Name of the person to greet.",
    ),
    add_on: str = typer.Option(
        "What a lovely day!",  # default value
        "--add_on",  # one way to call the option
        "--add-on",  # another way to call the option (what it would default to)
        help="Something else to say",  #
    ),
):
    """
    Greet a person with a simple hello message.

    Takes keyword arguments (typer.Option types), so order doesn't matter.

    Examples for CLI usage:
        >>> warp hello_world greet_keywords --name Andrew
        Hello, Andrew
        What a lovely day!
        >>> warp hello_world greet_keywords --name Andrew --add_on "Clever non sequitur"
        Hello, Andrew
        Clever non sequitur
        >>> warp hello_world greet_keywords --add_on "Clever non sequitur" --name Andrew
        Hello, Andrew
        Clever non sequitur

    """
    typer.echo(f"Hello, {name}!\n{add_on}")


@app.command(
    cls=CustomCLICommand,
    name="goodbye",
    short_help="Bid farewell to your comrade",
)
def goodbye(
    name: str = typer.Argument(..., help="Name of the person to bid farewell."),
):
    """
    Bid farewell to a person with a goodbye message.
    """
    typer.echo(f"Goodbye, {name}!")


@app.callback()
def callback():
    """
    Hello world module – courtesy of the imagewarp codebase.
    """
