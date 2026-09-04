#!/usr/bin/env python3
"""Reject retired source-project name tokens outside attribution surfaces."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

if not __debug__:
    raise RuntimeError(
        "this assertion-based check must not run with Python optimization"
    )

FORBIDDEN = ("pa" + "ri", "hec" + "ke")
RETIRED_NATIVE_IDENTITY = "silex" + "-" + "cpp"
PROJECT_STEM = "silex"
CPP_STEM = "cpp"
RETIRED_IDENTITY_PATTERNS = (
    (
        RETIRED_NATIVE_IDENTITY,
        re.compile(
            rf"(?i)(?:{PROJECT_STEM}[-_]{CPP_STEM}|"
            rf"{PROJECT_STEM}{CPP_STEM}|{PROJECT_STEM}\s+c\+\+|"
            rf"{PROJECT_STEM}\s*::\s*{CPP_STEM})"
        ),
    ),
)
UPSTREAM_CONTENT_PREFIX_ALLOWLIST = ("doc/", "docs/", "LICENSES/")
LEGAL_BASENAMES = {
    "COPYING",
    "COPYING.md",
    "COPYING.rst",
    "COPYING.txt",
    "LICENSE",
    "LICENSE.md",
    "LICENSE.rst",
    "LICENSE.txt",
    "NOTICE",
    "NOTICE.md",
    "NOTICE.rst",
    "NOTICE.txt",
    "THIRD_PARTY_NOTICES.md",
}
RETIRED_IDENTITY_CONTENT_ALLOWLIST = {"THIRD_PARTY_NOTICES.md"}


def word_components(text: str) -> list[str]:
    separated = re.sub(r"[^A-Za-z0-9]+", " ", text)
    separated = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", separated)
    separated = re.sub(r"(?<=[A-Z])(?=[A-Z][a-z])", " ", separated)
    separated = re.sub(r"(?<=[A-Za-z])(?=[0-9])", " ", separated)
    separated = re.sub(r"(?<=[0-9])(?=[A-Za-z])", " ", separated)
    return [part.casefold() for part in separated.split()]


def forbidden_components(text: str) -> set[str]:
    hits = set(word_components(text)).intersection(FORBIDDEN)
    runs = re.findall(r"[A-Za-z0-9]+", text)
    for run in runs:
        folded_run = run.casefold()
        for forbidden in FORBIDDEN:
            start = folded_run.find(forbidden)
            while start >= 0:
                end = start + len(forbidden)
                before = run[:start]
                after = run[end:]
                starts_run = start == 0
                follows_acronym = len(before) >= 2 and before.isupper()
                followed_by_digit = bool(after) and after[0].isdigit()
                followed_by_capitalized_word = (
                    len(after) >= 2 and after[0].isupper() and after[1].islower()
                )
                followed_by_numbered_acronym = (
                    bool(after)
                    and all(character.isupper() or character.isdigit() for character in after)
                    and any(character.isdigit() for character in after)
                )
                followed_by_acronym = len(after) >= 2 and after.isupper()
                if (
                    starts_run
                    and (
                        not after
                        or followed_by_digit
                        or followed_by_capitalized_word
                        or followed_by_acronym
                    )
                ) or (
                    follows_acronym
                    and (
                        followed_by_digit
                        or followed_by_capitalized_word
                        or followed_by_numbered_acronym
                    )
                ):
                    hits.add(forbidden)
                start = folded_run.find(forbidden, start + 1)
    for label, pattern in RETIRED_IDENTITY_PATTERNS:
        if pattern.search(text):
            hits.add(label)
    return hits


def allowed_content_hits(path: str) -> set[str]:
    name = Path(path).name
    allowed: set[str] = set()
    if (
        path == "AGENTS.md"
        or name == "README.md"
        or path.startswith(UPSTREAM_CONTENT_PREFIX_ALLOWLIST)
        or name in LEGAL_BASENAMES
    ):
        allowed.update(FORBIDDEN)
    if path in RETIRED_IDENTITY_CONTENT_ALLOWLIST:
        allowed.add(RETIRED_NATIVE_IDENTITY)
    return allowed


def tracked_worktree_paths(root: Path) -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    )
    return sorted(
        path.decode("utf-8")
        for path in completed.stdout.split(b"\0")
        if path
    )


def self_test() -> None:
    assert forbidden_components("pa" + "ri") == {FORBIDDEN[0]}
    assert forbidden_components("HEC" + "KE") == {FORBIDDEN[1]}
    assert forbidden_components("pa" + "ri_route") == {FORBIDDEN[0]}
    assert forbidden_components("Hec" + "keRoute") == {FORBIDDEN[1]}
    assert forbidden_components("XMLHec" + "keRoute") == {FORBIDDEN[1]}
    assert forbidden_components("XML" + FORBIDDEN[0].upper() + "2Route") == {
        FORBIDDEN[0]
    }
    assert forbidden_components(FORBIDDEN[0].upper() + "HTTP2") == {
        FORBIDDEN[0]
    }
    mixed_first = (
        FORBIDDEN[0][0]
        + FORBIDDEN[0][1].upper()
        + FORBIDDEN[0][2]
        + FORBIDDEN[0][3].upper()
    )
    mixed_second = (
        FORBIDDEN[1][0]
        + FORBIDDEN[1][1].upper()
        + FORBIDDEN[1][2]
        + FORBIDDEN[1][3].upper()
        + FORBIDDEN[1][4]
    )
    assert forbidden_components(mixed_first) == {FORBIDDEN[0]}
    assert forbidden_components(mixed_second) == {FORBIDDEN[1]}
    assert forbidden_components("src/" + mixed_first + "-route.cpp") == {
        FORBIDDEN[0]
    }
    assert forbidden_components("XML" + mixed_first + "2Route") == {
        FORBIDDEN[0]
    }
    assert forbidden_components("src/pa" + "ri-route.cpp") == {FORBIDDEN[0]}
    assert forbidden_components("Silex" + "Cpp") == {RETIRED_NATIVE_IDENTITY}
    assert forbidden_components("silex" + "_cpp_target") == {
        RETIRED_NATIVE_IDENTITY
    }
    assert forbidden_components("silex" + "-cpp-tool") == {
        RETIRED_NATIVE_IDENTITY
    }
    assert forbidden_components("Silex" + " C++") == {RETIRED_NATIVE_IDENTITY}
    assert forbidden_components("Silex::" + "cpp") == {RETIRED_NATIVE_IDENTITY}
    assert forbidden_components("Silex" + "CppConfig") == {
        RETIRED_NATIVE_IDENTITY
    }
    assert forbidden_components("lib" + "silex" + "_cpp.so") == {
        RETIRED_NATIVE_IDENTITY
    }
    assert not forbidden_components("Silex::silex")
    assert not forbidden_components("libsilex.so")
    assert not forbidden_components("comparison")
    assert not forbidden_components("inheritance")
    assert FORBIDDEN[0] in allowed_content_hits("LICENSE")
    assert FORBIDDEN[1] in allowed_content_hits("NOTICE.txt")
    assert FORBIDDEN[0] in allowed_content_hits("LICENSES/dependency.txt")
    assert RETIRED_NATIVE_IDENTITY not in allowed_content_hits("README.md")
    assert RETIRED_NATIVE_IDENTITY in allowed_content_hits(
        "THIRD_PARTY_NOTICES.md"
    )
    assert not allowed_content_hits("LICENSE_" + FORBIDDEN[0] + ".cpp")
    assert not allowed_content_hits("NOTICE_helper.py")

    with tempfile.TemporaryDirectory(prefix="source-name-boundary-") as temporary:
        root = Path(temporary)
        retired_path = "doc/" + RETIRED_NATIVE_IDENTITY + ".rst"
        ordinary_policy = "docs/development/testing.rst"
        upstream_policy = "docs/development/source_fidelity.rst"
        provenance = "THIRD_PARTY_NOTICES.md"
        fixtures = {
            retired_path: "neutral\n",
            ordinary_policy: "retired " + RETIRED_NATIVE_IDENTITY + " identity\n",
            "README.md": "retired " + RETIRED_NATIVE_IDENTITY + " identity\n",
            upstream_policy: "upstream " + FORBIDDEN[0] + " and " + FORBIDDEN[1] + "\n",
            provenance: "Imported from " + RETIRED_NATIVE_IDENTITY + " with attribution\n",
        }
        for relative, content in fixtures.items():
            fixture = root / relative
            fixture.parent.mkdir(parents=True, exist_ok=True)
            fixture.write_text(content, encoding="utf-8")
        policy_findings = scan_paths(root, sorted(fixtures))
        assert any(
            (retired_path + ": path contains " + RETIRED_NATIVE_IDENTITY)
            in finding
            for finding in policy_findings
        )
        assert any(
            (ordinary_policy + ":1: contains " + RETIRED_NATIVE_IDENTITY)
            in finding
            for finding in policy_findings
        )
        assert any(
            ("README.md:1: contains " + RETIRED_NATIVE_IDENTITY) in finding
            for finding in policy_findings
        )
        assert not any(
            finding.startswith(upstream_policy + ":") for finding in policy_findings
        )
        assert not any(
            finding.startswith(provenance + ":") for finding in policy_findings
        )

        safe_target = root / "safe-target"
        safe_target.write_text("neutral\n", encoding="utf-8")
        named_target = root / ("data-" + FORBIDDEN[0])
        named_target.write_text("neutral\n", encoding="utf-8")
        links = {
            "neutral-dangling": "missing-" + FORBIDDEN[0] + "-target",
            "neutral-live": named_target.name,
            FORBIDDEN[0] + "-dangling-link": "missing-safe-target",
            FORBIDDEN[0] + "-live-link": safe_target.name,
        }
        for name, target in links.items():
            os.symlink(target, root / name)
        findings = scan_paths(root, sorted(links))
        assert any("neutral-dangling:1: contains" in finding for finding in findings)
        assert any("neutral-live:1: contains" in finding for finding in findings)
        assert any(
            (FORBIDDEN[0] + "-dangling-link: path contains") in finding
            for finding in findings
        )
        assert any(
            (FORBIDDEN[0] + "-live-link: path contains") in finding
            for finding in findings
        )


def scan_paths(root: Path, relatives: list[str]) -> list[str]:
    findings: list[str] = []
    for relative in relatives:
        path = root / relative
        path_hits = forbidden_components(relative)
        if path_hits:
            findings.append(
                f"{relative}: path contains {','.join(sorted(path_hits))}"
            )
        if path.is_symlink():
            text = os.readlink(path)
            target_hits = forbidden_components(text)
            if target_hits:
                findings.append(
                    f"{relative}:1: contains {','.join(sorted(target_hits))}"
                )
            continue
        else:
            if not path.exists() or not path.is_file():
                continue
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
        allowed_hits = allowed_content_hits(relative)
        for line_number, line in enumerate(text.splitlines(), start=1):
            hits = forbidden_components(line) - allowed_hits
            if hits:
                findings.append(
                    f"{relative}:{line_number}: contains {','.join(sorted(hits))}"
                )
    return findings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    self_test()

    findings = scan_paths(root, tracked_worktree_paths(root))

    if findings:
        print("retired source-project name tokens found:", file=sys.stderr)
        for finding in findings:
            print(f"  {finding}", file=sys.stderr)
        return 1
    print("source-name boundary: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
