#!/usr/bin/env python3

"""
Regenerate the per-machine tables in the disruption parameters reference doc
from the `physics.attributes` section of each machine's `config.toml`.

Existing sections are replaced and new ones added, with all sections
sorted alphabetically by title. Content before the first section is left untouched.
"""

import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MACHINE_DIR = ROOT / "disruption_py" / "machine"
DOC = ROOT / "docs" / "usage" / "physics_methods" / "disruption_parameters_reference.md"

# Optional display names per machine directory (default: dir name)
NAMES = {
    "cmod": "C-Mod",
    "d3d": "DIII-D",
    "east": "EAST",
    "generic": "Generic",
    "hbtep": "HBT-EP",
    "mast": "MAST",
}


def format_num(x) -> str:
    """
    Format a number compactly for display in the validity range column.

    Uses general format (``%g``) and then strips the zero padding and ``+`` sign
    that Python adds to exponents, so ``5e+06`` is rendered as ``5e6``.

    Parameters
    ----------
    x : float or int
        Number to format.

    Returns
    -------
    str
        Compact representation of `x`.

    Notes
    -----
    General format keeps six significant digits and only switches to exponent
    notation for exponents below -5 or at least 6, so ``300000.0`` stays
    ``300000`` while ``5000000.0`` becomes ``5e6``. Non-finite values are passed
    through as ``inf``, ``-inf`` and ``nan``.

    Examples
    --------
    >>> format_num(5000000.0)
    '5e6'
    >>> format_num(0.2e20)
    '2e19'
    >>> format_num(1e-5)
    '1e-5'
    >>> format_num(0.025)
    '0.025'
    >>> format_num(float("inf"))
    'inf'
    """
    mantissa, _, exponent = f"{x:g}".partition("e")
    return f"{mantissa}e{int(exponent)}" if exponent else mantissa


def make_section(title: str, attributes: dict) -> str:
    """
    Build a markdown section containing a table of parameters.

    Parameters
    ----------
    title : str
        Section title, written as a level-2 heading (``## title``).
    attributes : dict
        Mapping of parameter name to its metadata, as found under
        ``physics.attributes`` in a machine ``config.toml``. Each value is a
        mapping that may contain the keys ``description``, ``units`` and
        ``validity``; ``validity`` is a ``[min, max]`` list of numbers.

    Returns
    -------
    str
        Markdown for the heading followed by a table with the columns Parameter,
        Description, Units and Validity Range, terminated by a newline.

    Notes
    -----
    Rows are sorted by parameter name. A missing ``description``, ``units`` or
    ``validity`` is shown as ``-``. Pipe characters in descriptions are escaped
    with a backslash so that they do not split the table cell. Validity bounds
    are formatted with `format_num`.
    """
    lines = [
        f"## {title}",
        "",
        "| Parameter | Description | Units | Validity Range |",
        "|---|---|---|---|",
    ]
    for name, attr in sorted(attributes.items()):
        desc = attr.get("description", "-").replace("|", "\\|")
        units = attr.get("units", "-")
        validity_list = attr.get("validity")
        validity = (
            f"[{', '.join(map(format_num, validity_list))}]" if validity_list else "-"
        )
        lines.append(f"| {name} | {desc} | {units} | {validity} |")
    return "\n".join(lines) + "\n"


def split_sections(text: str) -> tuple[str, dict[str, str]]:
    """
    Split a markdown document into a header and its level-2 sections.

    Parameters
    ----------
    text : str
        Full markdown text of the document.

    Returns
    -------
    header : str
        Everything before the first level-2 heading, unchanged.
    sections : dict of str to str
        Mapping of section title to the section text, including its heading line
        and trailing blank lines. Keys are in document order.

    Notes
    -----
    Only lines starting with ``"## "`` begin a new section, so deeper headings
    such as ``###`` stay inside the section that contains them. If a title
    appears more than once, only the text of its last occurrence is kept.
    """
    header, sections, title = [], {}, None
    for line in text.splitlines(keepends=True):
        if line.startswith("## "):
            title = line[3:].strip()
            sections[title] = ""
        if title is None:
            header.append(line)
        else:
            sections[title] += line
    return "".join(header), sections


def main():
    """
    Update the parameters reference doc from the machine ``config.toml`` files.

    Every subdirectory of ``disruption_py/machine`` that has a ``config.toml``
    is read. For each top-level table in the file that defines
    ``physics.attributes``, a section is generated and written into the doc,
    replacing the existing section with the same title or adding a new one. All
    sections are then sorted alphabetically by title, ignoring case, and the
    file is rewritten with the header preserved.

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        If the reference doc does not exist.
    tomllib.TOMLDecodeError
        If a ``config.toml`` is not valid TOML.

    Notes
    -----
    Section titles have the form ``"<name> Disruption Parameter Descriptions"``,
    where ``<name>`` comes from `NAMES`. A machine directory missing from
    `NAMES` falls back to its directory name and a note is printed. Machines
    without a ``config.toml`` or without ``physics.attributes`` are skipped.
    Sections for machines that no longer exist are not removed from the doc.
    """
    header, sections = split_sections(DOC.read_text(encoding="utf-8"))

    for machine_dir in sorted(MACHINE_DIR.iterdir()):
        config_file = machine_dir / "config.toml"
        if not config_file.is_file():
            continue
        with open(config_file, "rb") as f:
            toml = tomllib.load(f)
        for cfg in toml.values():
            attributes = cfg.get("physics", {}).get("attributes")
            if not attributes:
                continue
            name = NAMES.get(machine_dir.name)
            if name is None:
                name = machine_dir.name
                print(f"Note: '{machine_dir.name}' not in NAMES, using '{name}'")
            title = f"{name} Disruption Parameter Descriptions"
            sections[title] = make_section(title, attributes)
            print(f"{machine_dir.name}: {len(attributes)} parameters")

    # Case insensitive sorting
    ordered = sorted(sections.items(), key=lambda item: item[0].lower())
    body = "\n".join(text.strip("\n") + "\n" for _, text in ordered)
    DOC.write_text(header + body, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
