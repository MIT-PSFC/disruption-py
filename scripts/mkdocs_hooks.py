#!/usr/bin/env python3

"""
MkDocs hooks for the DisruptionPy documentation.

README.md and INSTALL.md link to each other with repository-relative paths,
which work on GitHub. The documentation pastes sections of those files into
pages under docs/ through pymdownx.snippets, where the same paths no longer
point anywhere. This hook rewrites such links to the documentation page that
carries the same content, after snippets are expanded and before MkDocs
resolves and validates links.

Example, for a section of README.md pasted into docs/installation/install.md:

    [Installation guide](INSTALL.md)
    -> [Installation guide](/installation/install.md)

The hook is registered in mkdocs.yml under `hooks:`, and MkDocs calls
`on_config` below once, when it loads the configuration.
"""

import re

from markdown import Extension
from markdown.preprocessors import Preprocessor

# Links to these repository files are rewritten.
# Keys: file names exactly as they appear in links in the repository's Markdown.
# Values: the documentation page with the same content, relative to docs/.
# To handle another file, add a line here; the pattern below picks it up.
REPO_FILE_TO_PAGE = {
    "INSTALL.md": "installation/install.md",
    "README.md": "index.md",
    "REFERENCES.md": "references.md",
}

# Matches the target part of a Markdown link to one of the files above,
# e.g. the `](INSTALL.md)` in `[Installation guide](INSTALL.md)`,
# or the `](INSTALL.md#public-installations)` in a link to a section.
# The file name must be the whole target, so links such as `](docs/README.md)`
# or `](https://.../README.md)` do not match.
# `INSTALL\.md|README\.md|...`: any one of the file names, dots escaped to match
# literally (an unescaped `.` would match any character)
_FILES = "|".join(re.escape(name) for name in REPO_FILE_TO_PAGE)
_LINK = re.compile(
    rf"""
    \]\(                       # `](`: end of the link text, start of the target
    (?P<file>{_FILES})         # one of the file names, e.g. `INSTALL.md`
    (?P<anchor>\#[^)\s]*)?     # optional section anchor, e.g. `#public-installations`
    \)                         # `)`: end of the target
    """,
    re.VERBOSE,  # allows the whitespace and comments above inside the pattern
)


def _rewrite_link(match: re.Match) -> str:
    """
    Return the replacement for one matched link target.

    The new target starts with `/`, i.e. it is relative to docs/: mkdocs.yml sets
    `validation: links: absolute_links: relative_to_docs`, so MkDocs then checks
    that the page exists and converts the target into a relative URL that works
    from whichever page the link ends up on.
    """
    page = REPO_FILE_TO_PAGE[match["file"]]
    anchor = match["anchor"] or ""
    return f"](/{page}{anchor})"


class RepoLinkPreprocessor(Preprocessor):
    """
    Rewrite links to repository files into links to documentation pages.

    A Markdown preprocessor receives a page's source as a list of lines, before
    it is converted to HTML, and returns the modified lines.
    """

    def run(self, lines: list[str]) -> list[str]:
        """Rewrite every matching link target on every line."""
        return [_LINK.sub(_rewrite_link, line) for line in lines]


class RepoLinkExtension(Extension):
    """Markdown extension that adds `RepoLinkPreprocessor` to the converter."""

    def extendMarkdown(self, md):  # pylint: disable=invalid-name
        """Register the preprocessor (called by Python-Markdown; name is fixed)."""
        # Preprocessors run from highest to lowest priority.
        # pymdownx.snippets, which pastes in the README/INSTALL sections, uses 32;
        # 31 makes this run right after it, so the pasted links are rewritten.
        # MkDocs resolves and validates links later, in a tree processor that
        # runs after all preprocessors.
        md.preprocessors.register(RepoLinkPreprocessor(md), "repo_links", 31)


def on_config(config):
    """
    MkDocs hook: add the extension to the configured Markdown extensions.

    MkDocs builds each page's Markdown converter from this list, so the
    preprocessor runs on every page.
    """
    config.markdown_extensions.append(RepoLinkExtension())
    return config
