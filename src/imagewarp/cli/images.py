import sys
from pathlib import Path
import subprocess
import re
from typing import Optional

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


# @app.command(
#     cls=CustomCLICommand,
#     name="run-tests",
#     short_help="Run pytest (optionally a single test) and emit its image.",
# )
# def run_tests(
#     nodeid: str = typer.Argument(
#         None,
#         help="(optional) pytest node id to run, e.g. "
#         "tests/test_utilities.py::test_checkerboard",
#     ),
#     quiet: bool = typer.Option(
#         False,
#         "--quiet",
#         "-q",
#         help="Suppress pytest-banner; emit only the image link.",
#     ),
# ):
#     """
#     Run pytest on either:

#       • the whole `tests/` tree (if no nodeid given), or
#       • exactly the nodeid you passed.

#     Any test wrapped in @cmdrun_test that writes into
#     docs/userguide/src/images will emit a Markdown image
#     link that cmdrun will pick up and render in your book.
#     """
#     # # Ensure our image folder exists #TODO hardcode fix
#     image_dir = Path("docs") / "userguide" / "src" / "images"
#     image_dir.mkdir(parents=True, exist_ok=True)

#     # Build pytest args
#     pytest_cmd = [
#         sys.executable,
#         "-m",
#         "pytest",
#         "-q",  # quiet pytest header
#         "-s",  # turn off capture so our print outs appear
#         "--tb=short",  # shorter tracebacks
#     ]

#     # If user passed a nodeid, run only that; otherwise run the whole suite
#     if nodeid:
#         pytest_cmd.append(nodeid)
#     else:
#         pytest_cmd.append("tests")

#     # Run from your project root so pytest can find tests/ and conftest.py
#     project_root = Path(__file__).parents[3]  # TODO hardcode fix?
#     book_src = project_root / "docs" / "userguide" / "src"  # TODO hardcode fix

#     result = subprocess.run(
#         pytest_cmd,
#         cwd=project_root,
#         capture_output=True,
#         text=True,
#     )

#     # optionally display pass/fail
#     if not quiet:
#         if result.returncode == 0:
#             typer.secho("✅ pytest succeeded", fg=typer.colors.GREEN)
#         else:
#             typer.secho("❌ pytest failed", fg=typer.colors.RED)

#     if result.returncode == 0:
#         # pull out only the markdown‐image line and rebase its path
#         # mdbook wants relative path, not absolute
#         for line in result.stdout.splitlines():
#             if line.startswith("!["):
#                 # extract the absolute path inside the parentheses
#                 m = re.match(r"!\[(.*?)\]\((/.*)\)", line)
#                 if m:
#                     alttext, abs_path = m.groups()
#                     abs_path = Path(abs_path)
#                     # turn it into a relative path under book_src:
#                     rel = abs_path.relative_to(book_src)
#                     # now echo a link mdbook can actually render:
#                     typer.echo(f"![{alttext}](images/{rel.name})")  # TODO hardcode
#                 else:
#                     # if it didn’t match, just echo it raw
#                     typer.echo(line)
#     else:
#         typer.secho(
#             "❌ pytest failed, image will not render", fg=typer.colors.RED
#         )  # TODO hardcode

#     # stderr in case pytest warnings/errors appeared
#     if result.stderr:
#         typer.echo(result.stderr, err=True)

#     raise typer.Exit(code=result.returncode)


@app.command(
    cls=CustomCLICommand,
    name="run-tests",
    short_help="Run pytest (optionally a single test) and emit a <figure>.",
)
def run_tests(
    nodeid: str = typer.Argument(
        ..., help="pytest node id, e.g. tests/test_utilities.py::test_checkerboard"
    ),
    caption: Optional[str] = typer.Option(
        None,
        "--caption",
        "-c",
        help="Caption for the figure (defaults to the alt-text).",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress pytest pass/fail banner; only emit <figure>.",
    ),
):
    """
    Run pytest on a single node. Any test wrapped in @cmdrun_test that
    writes into docs/userguide/src/images will have its image wrapped
    in a <figure>…</figure> block with an <img> and <figcaption>.
    """
    # ensure image folder exists
    project_root = Path(__file__).parents[3]  # TODO hardcode fix?
    book_src = project_root / "docs" / "userguide" / "src"  # TODO hardcode fix
    image_dir = book_src / "images"  # TODO hardcode fix
    image_dir.mkdir(parents=True, exist_ok=True)

    # build pytest invocation
    pytest_cmd = [sys.executable, "-m", "pytest", "-q", "-s", "--tb=short", nodeid]

    result = subprocess.run(
        pytest_cmd,
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    # optionally show pass/fail banner
    if not quiet:
        if result.returncode == 0:
            typer.secho("✅ pytest succeeded", fg=typer.colors.GREEN)
        else:
            typer.secho("❌ pytest failed", fg=typer.colors.RED)

    # now produce HTML figure(s)
    if result.returncode == 0:
        # scan stdout for lines like: ![alt](absolute/path/to/png)
        figures_emitted = False
        for line in result.stdout.splitlines():
            m = re.match(r"!\[(.*?)\]\((/.*?\.png)\)", line)
            if not m:
                continue

            alttext, abs_path = m.groups()
            abs_path = Path(abs_path)
            rel_path = abs_path.relative_to(book_src)

            # decide on figcaption
            figcap = caption if caption else alttext

            # emit the <figure> block
            typer.echo(
                "<figure>\n"
                f'  <img src="images/{rel_path.name}" alt="{alttext}">\n'
                f"  <figcaption>{figcap}</figcaption>\n"
                "</figure>"
            )
            figures_emitted = True

        if not figures_emitted:
            # no image line found
            typer.echo(
                "<figure>\n"
                f"  <figcaption>⚠️ No image generated by <code>{nodeid}</code>.</figcaption>\n"
                "</figure>"
            )

    else:
        # pytest failed
        typer.echo(
            "<figure>\n"
            f"  <figcaption>❌ pytest failed for <code>{nodeid}</code></figcaption>\n"
            "</figure>",
            err=False,
        )

    # always print stderr (so you see tracebacks if you inspect logs)
    if result.stderr:
        typer.echo(result.stderr, err=True)

    raise typer.Exit(code=result.returncode)


@app.callback()
def callback():
    """
    Image module – courtesy of the imagewarp codebase.
    """
