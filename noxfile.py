import nox
from laminci.nox import build_docs, run_pre_commit, run_pytest

nox.options.default_venv_backend = "none"


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
def build(session):
    session.run(*"pip install .[dev]".split())
    run_pytest(session)
    # lamin_sphinx imports lamindb while generating API docs. lamin-docs
    # installs lamindb and creates an instance before running lndocs.
    session.run(*"uv pip install --system lamindb".split())
    session.run(*"lamin init --storage ./docsbuild".split())
    build_docs(session)
