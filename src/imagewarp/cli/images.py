import sys
from pathlib import Path
import subprocess
import re

# import click  # unused import
import typer
from imagewarp.cli.customization import CustomCLIGroup, CustomCLICommand

app = typer.Typer(
    cls=CustomCLIGroup,
    no_args_is_help=True,
    help="Run tests for userguide image generation",
)

# Only force help at the module level when no subcommand is provided.
if len(sys.argv) == 2:
    sys.argv.append("--help")


@app.command(
    cls=CustomCLICommand,
    name="run-tests",
    short_help="Run pytest (optionally a single test) and emit its image.",
)
def run_tests(
    nodeid: str = typer.Argument(
        None,
        help="(optional) pytest node id to run, e.g. "
        "tests/test_utilities.py::test_checkerboard",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress pytest-banner; emit only the image link.",
    ),
):
    """
    Run pytest on either:

      • the whole `tests/` tree (if no nodeid given), or
      • exactly the nodeid you passed.

    Any test wrapped in @cmdrun_test that writes into
    docs/userguide/src/images will emit a Markdown image
    link that cmdrun will pick up and render in your book.
    """
    # Ensure our image folder exists
    image_dir = Path("docs") / "userguide" / "src" / "images"
    image_dir.mkdir(parents=True, exist_ok=True)

    # Build pytest args
    pytest_cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-q",  # quiet pytest header
        "-s",  # turn off capture so our print outs appear
        "--tb=short",  # shorter tracebacks
    ]

    # If user passed a nodeid, run only that; otherwise run the whole suite
    if nodeid:
        pytest_cmd.append(nodeid)
    else:
        pytest_cmd.append("tests")

    # Run from your project root so pytest can find tests/ and conftest.py
    project_root = Path(__file__).parents[3]
    book_src = project_root / "docs" / "userguide" / "src"

    result = subprocess.run(
        pytest_cmd,
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    # optionally display pass/fail
    if not quiet:
        if result.returncode == 0:
            typer.secho("✅ pytest succeeded", fg=typer.colors.GREEN)
        else:
            typer.secho("❌ pytest failed", fg=typer.colors.RED)

    # pull out only the markdown‐image line and rebase its path
    # mdbook wants relative path, not absolute
    for line in result.stdout.splitlines():
        if line.startswith("!["):
            # extract the absolute path inside the parentheses
            m = re.match(r"!\[(.*?)\]\((/.*)\)", line)
            if m:
                alttext, abs_path = m.groups()
                abs_path = Path(abs_path)
                # turn it into a relative path under book_src:
                rel = abs_path.relative_to(book_src)
                # now echo a link mdbook can actually render:
                typer.echo(f"![{alttext}](images/{rel.name})")
            else:
                # if it didn’t match, just echo it raw
                typer.echo(line)

    # stderr in case pytest warnings/errors appeared
    if result.stderr:
        typer.echo(result.stderr, err=True)

    raise typer.Exit(code=result.returncode)


# @app.command(
#     cls=CustomCLICommand,
#     name="run-tests",
#     short_help="Run unit tests",
# )
# def run_tests():
#     """
#     Run pytest on your utilities tests. Any images that the tests generate
#     into docs/userguide/src/images will be picked up by the cmdrun decorator
#     and printed out as Markdown image links.
#     """
#     # 1) Ensure the image directory exists
#     image_dir = Path("docs") / "userguide" / "src" / "images"
#     image_dir.mkdir(parents=True, exist_ok=True)

#     # 2) Build the pytest command
#     #    -m pytest ensures we use the same interpreter environment
#     pytest_args = [
#         sys.executable,
#         "-m",
#         "pytest",
#         "-q",  # quiet output
#         "--tb=short",  # shorter tracebacks
#         "tests/test_utilities.py",
#     ]

#     # 3) Run in the project root so pytest can discover tests and conftest
#     project_root = Path(__file__).parents[3]  # adjust if needed
#     result = subprocess.run(
#         pytest_args,
#         cwd=project_root,
#         capture_output=True,
#         text=True,
#     )

#     # 4) Echo results
#     if result.returncode == 0:
#         typer.secho("✅ All tests passed successfully.", fg=typer.colors.GREEN)
#     else:
#         typer.secho("❌ Some tests failed.", fg=typer.colors.RED)

#     # Always print stdout/stderr so cmdrun can pick up the decorator’s markdown
#     typer.echo(result.stdout)
#     typer.echo(result.stderr, err=True)

#     # Exit with the same code so CI knows if tests failed
#     raise typer.Exit(code=result.returncode)


@app.callback()
def callback():
    """
    Image module – courtesy of the imagewarp codebase.
    """
