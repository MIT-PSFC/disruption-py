#!/usr/bin/env python3

"""
Regenerate the per-machine tables in the disruption parameters reference doc
from the `physics.attributes` section of each machine's `config.toml`.

Existing sections are replaced in place, new ones are appended,
and any other content in the doc is left untouched.
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


def fmt_num(x) -> str:
    """Format a number compactly, e.g. 5000000.0 -> 5e6."""
    mantissa, _, exponent = f"{x:g}".partition("e")
    return f"{mantissa}e{int(exponent)}" if exponent else mantissa


def make_section(title: str, attributes: dict) -> str:
    """Build a markdown section with a table of parameters."""
    lines = [
        f"## {title}",
        "",
        "| Parameter | Description | Units | Validity Range |",
        "|---|---|---|---|",
    ]
    for name, attr in sorted(attributes.items()):
        desc = attr.get("description", "-").replace("|", "\\|")
        units = attr.get("units", "-")
        validity = attr.get("validity")
        valid = f"[{', '.join(map(fmt_num, validity))}]" if validity else "-"
        lines.append(f"| {name} | {desc} | {units} | {valid} |")
    return "\n".join(lines) + "\n"


def split_sections(text: str) -> tuple[str, dict[str, str]]:
    """Split a markdown doc into its preamble and `## ` sections keyed by title."""
    preamble, sections, title = [], {}, None
    for line in text.splitlines(keepends=True):
        if line.startswith("## "):
            title = line[3:].strip()
            sections[title] = ""
        if title is None:
            preamble.append(line)
        else:
            sections[title] += line
    return "".join(preamble), sections


def main():
    """Update the doc with one section per machine that defines attributes."""
    preamble, sections = split_sections(DOC.read_text(encoding="utf-8"))

    for machine_dir in sorted(MACHINE_DIR.iterdir()):
        config = machine_dir / "config.toml"
        if not config.is_file():
            continue
        with open(config, "rb") as f:
            data = tomllib.load(f)
        for cfg in data.values():
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

    body = "\n".join(s.strip("\n") + "\n" for s in sections.values())
    DOC.write_text(preamble + body, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
