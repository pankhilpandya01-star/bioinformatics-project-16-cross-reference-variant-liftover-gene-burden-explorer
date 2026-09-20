"""Audit the exact local publication candidate without creating a repository."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import uuid
from pathlib import Path
from urllib.parse import unquote

from .liftover import LiftoverError


EXCLUDED_PARTS = {
    ".git",
    ".pytest_cache",
    "__pycache__",
    "local",
    "build",
    "dist",
}
TEXT_SUFFIXES = {
    ".csv",
    ".md",
    ".py",
    ".toml",
    ".txt",
    ".yml",
    ".yaml",
    ".bed",
    ".paf",
    ".fasta",
    ".fai",
    ".vcf",
    ".json",
    ".gitattributes",
    ".gitignore",
}
MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def candidate_files(root: Path, audit_output: Path) -> list[Path]:
    files: list[Path] = []
    audit_output = audit_output.resolve()
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        resolved = path.resolve()
        if audit_output == resolved or audit_output in resolved.parents:
            continue
        relative = path.relative_to(root)
        if any(part in EXCLUDED_PARTS or part.endswith(".egg-info") for part in relative.parts):
            continue
        files.append(path)
    return sorted(files, key=lambda item: item.relative_to(root).as_posix())


def _read_public_text(path: Path) -> str | None:
    lower_name = path.name.lower()
    if lower_name.endswith((".vcf.gz", ".csv.gz")):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    if path.suffix.lower() in TEXT_SUFFIXES or path.name in {"LICENSE"}:
        return path.read_text(encoding="utf-8", errors="replace")
    return None


def audit_text(path: Path, content: str) -> list[str]:
    findings: list[str] = []
    lowered = content.lower()
    forbidden_literals = {
        "native Windows user path": "c:" + "\\users\\",
        "mounted Windows user path": "/mnt/" + "c/users/",
        "synchronized-folder path": "one" + "drive",
        "unwanted authorship term": "chat" + "gpt",
        "unwanted vendor term": "open" + "ai",
        "unwanted generation term": "ai" + "-generated",
        "unwanted authorship phrase": "artificial " + "intelligence",
        "private-key marker": "begin " + "private key",
    }
    for label, value in forbidden_literals.items():
        if value in lowered:
            findings.append(f"{label}: {path}")
    secret_assignment = re.compile(
        r"(?i)(?:api[_-]?key|password|client[_-]?secret|access[_-]?token)\s*[:=]\s*['\"][^'\"]+"
    )
    if secret_assignment.search(content):
        findings.append(f"credential-like assignment: {path}")
    if re.search(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", content):
        findings.append(f"email address: {path}")
    return findings


def broken_markdown_links(root: Path, files: list[Path]) -> list[str]:
    broken: list[str] = []
    for path in files:
        if path.suffix.lower() != ".md":
            continue
        content = path.read_text(encoding="utf-8", errors="replace")
        for raw_target in MARKDOWN_LINK.findall(content):
            target = raw_target.strip().strip("<>")
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            target = unquote(target.split("#", 1)[0])
            if not target:
                continue
            resolved = (path.parent / target).resolve()
            if not resolved.exists():
                broken.append(f"{path.relative_to(root).as_posix()} -> {raw_target}")
    return broken


def publish_audit(root: Path, output_dir: Path) -> dict[str, object]:
    root = root.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise LiftoverError("publication-audit output directory must not already exist")
    if root not in output_dir.parents:
        raise LiftoverError("publication-audit output must be inside the project root")
    stage = output_dir.parent / f".{output_dir.name}.staging-{uuid.uuid4().hex}"
    stage.mkdir(parents=True)
    try:
        files = candidate_files(root, output_dir)
        if not files:
            raise LiftoverError("publication candidate contains no files")
        inventory: list[dict[str, object]] = []
        findings: list[str] = []
        for path in files:
            relative = path.relative_to(root).as_posix()
            size = path.stat().st_size
            inventory.append(
                {
                    "path": relative,
                    "bytes": size,
                    "sha256": sha256(path),
                }
            )
            content = _read_public_text(path)
            if content is not None:
                findings.extend(audit_text(Path(relative), content))
        broken = broken_markdown_links(root, files)
        if findings or broken:
            detail = "; ".join(findings + [f"broken link: {item}" for item in broken])
            raise LiftoverError(f"publication audit failed: {detail}")
        largest = max(inventory, key=lambda item: int(item["bytes"]))
        if int(largest["bytes"]) >= 100 * 1024 * 1024:
            raise LiftoverError(f"public file reaches GitHub's 100 MB limit: {largest['path']}")
        with (stage / "file_inventory.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["path", "bytes", "sha256"], lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(inventory)
        summary = {
            "workflow": "variant-liftover-explorer",
            "milestone": "publication_audit",
            "status": "pass",
            "candidate_files": len(inventory),
            "candidate_bytes": sum(int(item["bytes"]) for item in inventory),
            "largest_file": largest,
            "individual_file_limit_bytes": 100 * 1024 * 1024,
            "markdown_links_checked": True,
            "broken_markdown_links": 0,
            "privacy_and_credential_findings": 0,
            "authorship_wording_findings": 0,
            "local_storage_excluded": True,
            "repository_created": (root / ".git").is_dir(),
            "inventory_sha256": sha256(stage / "file_inventory.csv"),
        }
        (stage / "audit_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(stage, output_dir)
        return summary
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
