"""
This module provides a custom command-line interface (CLI) toolkit,
utilizing the Typer and Click libraries for enhanced user interaction
and command management.

Key Features:
- Custom CLI Group and Command classes that extend Typer and Click
functionalities to provide rich help formatting, including usage
instructions, parameter details, and error handling.
- A color-coded output system using the Rich library to improve the
readability of command-line messages and help texts.
- A logo display for branding purposes when the main command is
invoked.
- Enhanced error handling that provides user-friendly messages and
suggestions for available commands when invalid inputs are given.

Classes:
- Color: An enumeration of color styles used for formatting CLI output.
- CustomCLIGroup: A custom group class that manages subcommands,
formats help messages, and handles command resolution.
- CustomCLICommand: A custom command class that manages individual
command behavior, including parameter parsing and help display.

Constants:
- LOGO_TEXT: A string constant that contains the ASCII art logo for the
application.

Usage:
To use this module, instantiate the CustomCLIGroup and CustomCLICommand
classes to define your CLI structure, and utilize the provided color
formatting and help display features to enhance user experience.

Example:
```python
app = typer.Typer(cls=CustomCLIGroup)
@app.command(cls=CustomCLICommand)
def example_command(param: str):
    '''
    An example command that demonstrates the usage of the custom CLI.
    '''
    pass

"""

import sys
from enum import Enum
from enum import StrEnum
from typing import Final

import click
import typer
from rich.console import Console
from rich.table import Table
from rich import box

LOGO_TEXT: Final = r"""
    ██╗███╗   ███╗ █████╗  ██████╗ ███████╗ ░░╗     ░░╗ ░░░░░╗ ░░░░░░╗ ░░░░░░╗ 
    ██║████╗ ████║██╔══██╗██╔════╝ ██╔════╝ ░░║ ░░  ░░║░░╔══░░╗░░╔══░░╗░░╔══░░╗
    ██║██╔████╔██║███████║██║  ███╗█████╗   ░░╚░░░░╔░░║░░░░░░░║░░░░░░╔╝░░░░░░╔╝
    ██║██║╚██╔╝██║██╔══██║██║   ██║██╔══╝   ░░░░╝ ░░░░║░░╔══░░║░░╔══░░╗░░╔═══╝ 
    ██║██║ ╚═╝ ██║██║  ██║╚██████╔╝███████╗ ░░░╝   ░░░╗░░║  ░░║░░║  ░░║░░║ 
    ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝

                             ~~ ImageWarp ~~
                       Transform • Distort • Play
"""

console = Console()


class Color(StrEnum):
    """Colors for custom typer formmatting

    full list of colors for customization here:
    https://rich.readthedocs.io/en/stable/appendix/colors.html
    """

    FRONT_MATTER = "dark_sea_green1"
    MODULE = "magenta"
    COMMAND = "orange1"
    ERROR = "red"
    STATEMENT = "white"
    CURRENT_LEVEL = "yellow"
    HEADER = "bold cyan"
    PARAMETER = "green"
    TYPE = "blue"
    REQUIRED = "yellow"
    DEAFULT = "magenta"
    SHORT_HELP = "white"


class CustomCLIGroup(typer.core.TyperGroup):
    """
    CustomCLIGroup is an enhanced command group for module-level help
    and error handling in the CLI.

    This class extends Typer's base group functionality to provide a
    rich, color-coded help and usage display using the Rich library. It
    is designed to improve user interaction by offering:

      1. A dynamic "Usage:" line that adapts based on the current
         command context.
      2. Display of a short help message (if available) for the current
         command group.
      3. A formatted table listing subcommands (or modules) with
         color-coded entries:
           - Module names are shown in a designated module color.
           - Commands are shown in a designated command color.
      4. Full help documentation (docstring) display when available.

    Additional Features:
      - When the main command is invoked without additional subcommands,
        an ASCII art logo is printed to brand the CLI.
      - The resolve_command method is overridden to catch cases where a
        user supplies a non-existent command or module. In such cases,
        it prints a descriptive error message, indicates available
        options, and then displays the group help.
      - Overrides several methods (such as format_help, get_help, main)
        to integrate rich formatting and improved error management.

    Methods:
      list_commands(ctx):
          Returns a sorted list of available command names.
      format_help(ctx, formatter):
          Displays the help message for the current group.
      get_help(ctx):
          Displays help and exits the program.
      show_help(ctx):
          Formats and prints the usage, short help, subcommand table,
          and full help for the group.
      format_exception(ctx, exc):
          Overrides exception formatting (returns an empty string, as
          errors are handled separately).
      resolve_command(ctx, args):
          Attempts to resolve a command; on failure, prints an error
          message, shows available options, and exits.
      main(*args, **kwargs):
          Executes the main command loop with standalone mode disabled
          to allow error bubbling.
    """

    def list_commands(self, _ctx):
        """List all the commands (functions) added to the module"""
        return sorted(self.commands.keys())

    def format_help(self, ctx, _formatter):
        """Displays the help message for the current group."""
        self.show_help(ctx)

    def get_help(self, ctx):
        """Displays help and exits the program to handle errors."""
        self.show_help(ctx)
        ctx.exit()
        return ""

    def show_help(self, ctx: click.Context):
        """Custom show help function displaying table for groups"""
        # 1. Print usage header.
        # Determine the usage line based on the current command path:
        if ctx.command_path.strip() == "warp":
            console.print(f"[{Color.FRONT_MATTER}]{LOGO_TEXT}[/{Color.FRONT_MATTER}]")
            # If only "warp" is provided, show [MODULE]:
            usage_line = (
                f"[{Color.CURRENT_LEVEL}]warp[/{Color.CURRENT_LEVEL}] "
                f"[[{Color.MODULE}]MODULE[/{Color.MODULE}]] "
                f"[[{Color.COMMAND}]COMMAND[/{Color.COMMAND}]] "
                f"[[{Color.PARAMETER}]PARAMETERS[/{Color.PARAMETER}]]"
            )
            column_label = "Module"
        else:
            # Otherwise, just append [COMMAND] [OPTIONS]:
            usage_line = (
                f"[{Color.CURRENT_LEVEL}]{ctx.command_path}[/{Color.CURRENT_LEVEL}] "
                f"[[{Color.COMMAND}]COMMAND[/{Color.COMMAND}]] "
                f"[[{Color.PARAMETER}]PARAMETERS[/{Color.PARAMETER}]]"
            )
            column_label = "Command"
        console.print(f"[bold]Usage:[/] {usage_line}\n")

        # 2. Print short help for this group if available.
        short_help = getattr(self, "short_help", "")
        if short_help:
            console.print(short_help + "\n", style=f"{Color.SHORT_HELP}")
        # 3. Build a table listing subcommands.
        table = Table(show_header=True, header_style=f"{Color.HEADER}", box=box.MINIMAL)
        table.add_column(column_label, justify="left")
        table.add_column("Short help", justify="left")
        for name, cmd in self.commands.items():
            # Distinguish modules from commands:
            if isinstance(cmd, typer.core.TyperGroup):
                # For modules, display the name in {Color.MODULE}.
                display_name = f"[{Color.MODULE}]{name}[/{Color.MODULE}]"
            else:
                # For function commands, display the name in {Color.COMMAND}.
                display_name = f"[{Color.COMMAND}]{name}[/{Color.COMMAND}]"
            short = getattr(cmd, "short_help", "") or ""
            table.add_row(display_name, short)
        console.print(table)
        # 4. Print full help (docstring) if provided.
        if self.help:
            console.print("[bold]Full Help:[/]")
            console.print("\n" + self.help)

    def format_exception(self, _ctx, _exc):
        """Format exceptions generated by typer"""
        # We override this method but do not print anything here.
        return ""

    def resolve_command(self, ctx, args):
        """Attempts to resolve a command; on failure, prints an error
        message, shows available options, and exits."""
        try:
            return super().resolve_command(ctx, args)
        except click.exceptions.UsageError as e:
            # Print error message.
            console.print(f"Error: {e}\n", style=f"{Color.ERROR}")
            # Print a message indicating viable options.
            console.print(
                f"\n[{Color.STATEMENT}]Availble options at this level ([/{Color.STATEMENT}]"
                f"[{Color.CURRENT_LEVEL}]{ctx.command_path}[/{Color.CURRENT_LEVEL}]"
                f"[{Color.STATEMENT}]) are:\n[/{Color.STATEMENT}]"
            )
            # Show this group's help then exit.
            self.show_help(ctx)
            ctx.exit(1)
            return None  # This return is never reached but silences pylint.

    def main(self, *args, **kwargs):
        """Executes the main command loop with standalone mode disabled
        to allow error bubbling."""
        # Force standalone_mode=False so errors bubble up.
        kwargs.setdefault("standalone_mode", False)
        try:
            return super().main(*args, **kwargs)
        except Exception as e:  # pylint: disable=broad-exception-caught
            # catch broad exceptions
            console.print(f"Error: {e}\n", style=f"{Color.ERROR}")
            sys.exit(1)


class CustomCLICommand(click.Command):
    """
    CustomCLICommand is an enhanced command class for individual CLI
    commands that provides detailed, color-formatted help and error
    handling using the Rich library.

    This class extends Click's Command functionality to:
      1. Print a "Usage:" line with the current command path and its
         parameters.
      2. Display the command's short help message (if provided) in a
         highlighted style.
      3. Create a formatted table of parameters, detailing each
         parameter's name, type, whether it is required, its default
        value, and any associated short help.
      4. Output the full detailed help (the command's docstring) when
         available.

    In addition, CustomCLICommand intercepts missing parameter
    exceptions during argument parsing or command invocation. When such
    exceptions occur, it prints an informative error message, displays
    the help for the command, and exits the execution, ensuring that
    users are guided toward correct usage.

    Methods:
        parse_args(ctx, args):
            Parses the command-line arguments and, if a required
            parameter is missing, prints an error message and shows the
            command's help before exiting.
        invoke(ctx):
            Invokes the command; if a missing parameter is encountered
            during invocation, handles the error similarly to
            parse_args.
        format_help(ctx, formatter):
            Displays the help message for the command by calling the
            custom show_help method.
        get_help(ctx):
            Displays the help, then exits the process.
        show_help(ctx):
            Formats and prints:
              - A usage header
              - The command's short help (if available)
              - A detailed table of parameters
              - The full detailed help (if provided)
        format_exception(ctx, exc):
            Returns an empty string to allow higher-level error
            handling.
    """

    def __init__(self, *args, **kwargs):
        kwargs.pop("rich_markup_mode", None)
        kwargs.pop("rich_help_panel", None)
        super().__init__(*args, **kwargs)

    def parse_args(self, ctx, args):
        try:
            return super().parse_args(ctx, args)
        except click.MissingParameter as e:
            console.print(f"Error: {e}\n", style=f"{Color.ERROR}")
            self.show_help(ctx)
            ctx.exit(1)
            return None  # This return is never reached but silences pylint.

    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except click.MissingParameter as e:
            console.print(f"Error: {e}\n", style=f"{Color.ERROR}")
            self.show_help(ctx)
            ctx.exit(1)
            return None  # This return is never reached but silences pylint.

    def format_help(self, ctx, _formatter):
        self.show_help(ctx)

    def get_help(self, ctx):
        self.show_help(ctx)
        ctx.exit()
        return ""

    def show_help(self, ctx: click.Context):
        """Custom show help function displaying table for commands"""
        # 1. Print usage header for the command.
        console.print(
            f"[bold]Usage:[/] "
            f"[{Color.CURRENT_LEVEL}]{ctx.command_path}[/{Color.CURRENT_LEVEL}]  "
            f"[[{Color.PARAMETER}]PARAMETERS[/{Color.PARAMETER}]]\n"
        )

        # 2. Print short help if available.
        short_help = getattr(self, "short_help", "")
        if short_help:
            console.print(short_help + "\n", style=f"{Color.SHORT_HELP}")
        # 3. Build table for parameters with additional detail.
        table = Table(show_header=True, header_style=f"{Color.HEADER}", box=box.MINIMAL)
        table.add_column("Parameter", style=f"{Color.PARAMETER}", justify="left")
        table.add_column("Type", style=f"{Color.TYPE}", justify="left")
        table.add_column("Required", style=f"{Color.REQUIRED}", justify="center")
        table.add_column("Default", style=f"{Color.DEAFULT}", justify="left")
        table.add_column("Short help", style=f"{Color.SHORT_HELP}", justify="left")
        for param in self.params:
            if isinstance(param, click.Argument):
                param_type = "Argument"
                required = "Yes" if param.required else "No"
                default = str(param.default) if param.default is not None else ""
                help_text = getattr(param, "help", "")
                table.add_row(param.name, param_type, required, default, help_text)
            elif isinstance(param, click.Option):
                param_type = "Option"
                required = "Yes" if param.required else "No"
                default = str(param.default) if param.default is not None else ""
                help_text = param.help or ""
                names = ", ".join(param.opts)
                table.add_row(names, param_type, required, default, help_text)
        console.print("[bold]Parameters:[/]")
        if table.row_count:
            console.print(table)
        else:
            console.print("No parameters.\n")
        # 4. Print full detailed help (docstring) if provided.
        if self.help:
            console.print("[bold]Full Help:[/]")
            console.print("\n" + self.help)

    def format_exception(self, _ctx, _exc):
        """Format exceptions generated by typer"""
        # Do nothing here so that invoke/parse_args handles it.
        return ""
