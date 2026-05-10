#!/usr/bin/env python3

from __future__ import annotations

import argparse
import copy
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Any

try:
    import yaml
except ModuleNotFoundError as exc:
    raise SystemExit(
        "fatal: PyYAML is required to run integration suites. "
        "Install it with `python3 -m pip install pyyaml`."
    ) from exc


ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
TEMPLATE_RE = re.compile(r"\{([^{}]+)\}")
LINE_PATTERN_RE = re.compile(r"%\(([^()]+)\)")
DEFAULT_OVERLAY = {
    "defaults": {
        "env": {
            "LC_ALL": "C",
            "OMP_NUM_THREADS": "1",
            "MP_BINPATH": "{scratch}/tmp",
        }
    }
}
BIN_PATH_SUBDIRS = (
    "",
    "bin",
    "tools",
    "models",
    "bin/tools",
    "bin/models",
)
MULTI_CONFIG_SUBDIRS = (
    "Debug",
    "Release",
    "RelWithDebInfo",
    "MinSizeRel",
)


class SuiteError(RuntimeError):
    pass


def strip_ansi(text: str) -> str:
    return ANSI_ESCAPE_RE.sub("", text)


def is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def normalize_config(config: str | None) -> str | None:
    if config is None:
        return None
    config = config.strip()
    if not config:
        return None
    for known_config in MULTI_CONFIG_SUBDIRS:
        if config.lower() == known_config.lower():
            return known_config
    return config


def ordered_multi_config_subdirs(config: str | None) -> list[str]:
    config = normalize_config(config)
    if config is None:
        return list(MULTI_CONFIG_SUBDIRS)
    return [config] + [
        known_config for known_config in MULTI_CONFIG_SUBDIRS
        if known_config != config
    ]


def executable_search_dirs(root: Path, config: str | None = None) -> list[Path]:
    dirs: list[Path] = []
    seen: set[str] = set()
    active_config = normalize_config(config)
    config_subdirs = ordered_multi_config_subdirs(config)

    for subdir in BIN_PATH_SUBDIRS:
        base = root / subdir if subdir else root
        if active_config is None:
            candidates = [base]
            candidates.extend(base / config_subdir for config_subdir in config_subdirs)
        else:
            candidates = [base / config_subdir for config_subdir in config_subdirs]
            candidates.append(base)

        for candidate in candidates:
            key = str(candidate.resolve(strict=False))
            if key in seen:
                continue
            seen.add(key)
            dirs.append(candidate)

    return dirs


def deep_copy(value: Any) -> Any:
    return copy.deepcopy(value)


def resolve_reference(expr: str, context: dict[str, Any]) -> Any:
    if "." not in expr and "bindings" in context and expr in context["bindings"]:
        return context["bindings"][expr]
    current: Any = context
    for part in expr.split("."):
        if isinstance(current, dict) and part in current:
            current = current[part]
            continue
        raise SuiteError(f"Unknown template reference {{{expr}}}")
    return current


def render_value(value: Any, context: dict[str, Any]) -> Any:
    if isinstance(value, str):
        def replace(match: re.Match[str]) -> str:
            resolved = resolve_reference(match.group(1), context)
            return str(resolved)

        return TEMPLATE_RE.sub(replace, value)
    if isinstance(value, list):
        return [render_value(item, context) for item in value]
    if isinstance(value, dict):
        return {key: render_value(item, context) for key, item in value.items()}
    return value


def parse_float(text: str) -> float:
    return float(text.strip())


def parse_complex(text: str) -> complex:
    normalized = text.strip().replace(" ", "")
    if normalized.endswith("i"):
        normalized = normalized[:-1] + "j"
    normalized = normalized.replace("+j", "+1j").replace("-j", "-1j")
    if normalized == "j":
        normalized = "1j"
    if normalized == "-j":
        normalized = "-1j"
    return complex(normalized)


def extract_json_path(value: Any, path: str) -> Any:
    current = value
    if not path:
        return current
    for part in path.split("."):
        match = re.fullmatch(r"([^\[\]]+)(?:\[(\d+)\])?", part)
        if not match:
            raise SuiteError(f"Invalid json path segment: {part}")
        key = match.group(1)
        index_text = match.group(2)
        if not isinstance(current, dict) or key not in current:
            raise SuiteError(f"Missing json key while resolving path {path}: {key}")
        current = current[key]
        if index_text is not None:
            if not isinstance(current, list):
                raise SuiteError(f"Json path {path} expected list at {key}")
            index = int(index_text)
            try:
                current = current[index]
            except IndexError as exc:
                raise SuiteError(f"Json path {path} index out of range: {index}") from exc
    return current


def coerce_expected(value: Any, actual: Any) -> Any:
    if isinstance(actual, complex):
        if isinstance(value, (int, float)):
            return complex(value, 0.0)
        if isinstance(value, str):
            return parse_complex(value)
    if isinstance(actual, float):
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            return parse_float(value)
    return value


def approx_equal(actual: Any, expected: Any, abs_tol: float, rel_tol: float) -> bool:
    if isinstance(actual, complex) or isinstance(expected, complex):
        actual_complex = complex(actual)
        expected_complex = complex(expected)
        delta = abs(actual_complex - expected_complex)
        scale = max(abs(actual_complex), abs(expected_complex))
        return delta <= max(abs_tol, rel_tol * scale)
    return math.isclose(float(actual), float(expected), abs_tol=abs_tol, rel_tol=rel_tol)


def normalize_capture(text: str) -> str:
    return strip_ansi(text).replace("\r\n", "\n").replace("\r", "\n")


def line_pattern_literal_to_regex(literal: str) -> str:
    parts: list[str] = []
    i = 0
    while i < len(literal):
        if literal[i].isspace():
            while i < len(literal) and literal[i].isspace():
                i += 1
            parts.append(r"\s*")
            continue
        parts.append(re.escape(literal[i]))
        i += 1
    return "".join(parts)


def line_pattern_type_regex(kind: str) -> str:
    if kind == "float":
        return r"[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?"
    if kind == "int":
        return r"[+-]?\d+"
    if kind == "complex":
        return r"\S+"
    if kind == "word":
        return r"\S+"
    if kind == "rest":
        return r".+"
    raise SuiteError(f"Unsupported line-pattern placeholder type: {kind}")


def parse_line_pattern_value(kind: str, text: str) -> Any:
    if kind == "float":
        return parse_float(text)
    if kind == "int":
        return int(text)
    if kind == "complex":
        return parse_complex(text)
    if kind in {"word", "rest"}:
        return text
    raise SuiteError(f"Unsupported line-pattern placeholder type: {kind}")


def compile_line_pattern(pattern: str) -> tuple[re.Pattern[str], list[dict[str, Any]]]:
    pieces: list[str] = []
    captures: list[dict[str, Any]] = []
    position = 0

    for index, match in enumerate(LINE_PATTERN_RE.finditer(pattern)):
        pieces.append(line_pattern_literal_to_regex(pattern[position:match.start()]))

        spec = match.group(1).strip()
        if not spec:
            raise SuiteError(f"Invalid empty line-pattern placeholder in {pattern!r}")

        if ":" in spec:
            kind, name = spec.split(":", 1)
            kind = kind.strip()
            name = name.strip()
        else:
            kind = spec
            name = str(index)

        if not name:
            name = str(index)

        group_name = f"capture_{index}"
        pieces.append(f"(?P<{group_name}>{line_pattern_type_regex(kind)})")
        captures.append({"group": group_name, "kind": kind, "name": name})
        position = match.end()

    pieces.append(line_pattern_literal_to_regex(pattern[position:]))
    regex = re.compile(r"^" + "".join(pieces) + r"$")
    return regex, captures


def extract_line_pattern(text: str, pattern: str, capture: Any = None) -> Any:
    regex, captures = compile_line_pattern(pattern)

    for raw_line in normalize_capture(text).splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = regex.fullmatch(line)
        if match is None:
            continue

        ordered_values: list[Any] = []
        named_values: dict[str, Any] = {}
        for item in captures:
            value = parse_line_pattern_value(item["kind"], match.group(item["group"]))
            ordered_values.append(value)
            named_values[item["name"]] = value

        if capture is None:
            if not ordered_values:
                return line
            if len(ordered_values) == 1:
                return ordered_values[0]
            return named_values

        if isinstance(capture, int):
            try:
                return ordered_values[capture]
            except IndexError as exc:
                raise SuiteError(
                    f"Line-pattern capture index {capture} out of range for pattern {pattern!r}"
                ) from exc

        capture_text = str(capture)
        if capture_text.isdigit() or (
            capture_text.startswith("-") and capture_text[1:].isdigit()
        ):
            index = int(capture_text)
            try:
                return ordered_values[index]
            except IndexError as exc:
                raise SuiteError(
                    f"Line-pattern capture index {index} out of range for pattern {pattern!r}"
                ) from exc

        if capture_text not in named_values:
            raise SuiteError(
                f"Line-pattern capture {capture_text!r} not found in pattern {pattern!r}"
            )
        return named_values[capture_text]

    raise SuiteError(f"No output line matched pattern {pattern!r}")


def last_nonblank_line(text: str) -> str:
    lines = [line.strip() for line in normalize_capture(text).splitlines() if line.strip()]
    if not lines:
        raise SuiteError("Expected output, but command produced no non-blank lines")
    return lines[-1]


def indent_block(text: str, prefix: str = "  ") -> str:
    stripped = text.rstrip("\n")
    if not stripped:
        return prefix + "(empty)"
    return "\n".join(prefix + line for line in stripped.splitlines())


def snapshot_tree(root: Path, ignore_roots: list[Path] | None = None) -> set[Path]:
    root = root.resolve()
    ignored = [path.resolve() for path in (ignore_roots or []) if path.exists()]
    snapshot: set[Path] = set()

    for dirpath, dirnames, filenames in os.walk(root):
        current = Path(dirpath).resolve()
        dirnames[:] = [
            name
            for name in dirnames
            if not any(is_relative_to(current / name, ignore_root) for ignore_root in ignored)
        ]

        if current != root and not any(is_relative_to(current, ignore_root) for ignore_root in ignored):
            snapshot.add(current.relative_to(root))

        for filename in filenames:
            path = (current / filename).resolve()
            if any(is_relative_to(path, ignore_root) for ignore_root in ignored):
                continue
            snapshot.add(path.relative_to(root))

    return snapshot


def deep_merge(base: Any, overlay: Any) -> Any:
    if isinstance(base, dict) and isinstance(overlay, dict):
        merged = deep_copy(base)
        for key, value in overlay.items():
            if key in merged:
                merged[key] = deep_merge(merged[key], value)
            else:
                merged[key] = deep_copy(value)
        return merged
    return deep_copy(overlay)


def load_raw_suite(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    merged: dict[str, Any] = {}
    for import_name in data.get("imports", []):
        import_path = (path.parent / import_name).resolve()
        merged = deep_merge(merged, load_raw_suite(import_path))

    local_data = deep_copy(data)
    local_data.pop("imports", None)
    return deep_merge(merged, local_data)


COMPARE_KEYS = {"approx", "equals", "matches", "contains", "exists", "abs", "rel", "compare"}


def build_compare(spec: dict[str, Any]) -> dict[str, Any]:
    if "compare" in spec:
        return deep_copy(spec["compare"])
    if "approx" in spec:
        compare = {"kind": "approx", "value": deep_copy(spec["approx"])}
        if "abs" in spec:
            compare["abs"] = deep_copy(spec["abs"])
        if "rel" in spec:
            compare["rel"] = deep_copy(spec["rel"])
        return compare
    if "equals" in spec:
        return {"kind": "equals", "value": deep_copy(spec["equals"])}
    if "matches" in spec:
        return {"kind": "matches", "value": deep_copy(spec["matches"])}
    if "contains" in spec:
        return {"kind": "contains", "value": deep_copy(spec["contains"])}
    if "exists" in spec:
        return {"kind": "exists", "value": deep_copy(spec["exists"])}
    raise SuiteError("Assertion is missing a comparison")


def strip_compare_fields(spec: dict[str, Any]) -> dict[str, Any]:
    return {key: deep_copy(value) for key, value in spec.items() if key not in COMPARE_KEYS}


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def available_action_names(recipes: dict[str, Any]) -> set[str]:
    return {
        name
        for name, spec in recipes.items()
        if isinstance(spec, dict) and spec.get("kind") == "action"
    }


def available_probe_names(recipes: dict[str, Any]) -> set[str]:
    return {
        name
        for name, spec in recipes.items()
        if isinstance(spec, dict) and spec.get("kind") == "probe"
    }


def resolve_action_spec(op_name: str, recipes: dict[str, Any]) -> dict[str, Any]:
    if op_name in recipes and recipes[op_name].get("kind") == "action":
        return deep_copy(recipes[op_name])
    raise SuiteError(f"Unknown action operation: {op_name}")


def resolve_probe_spec(op_name: str, recipes: dict[str, Any]) -> dict[str, Any]:
    if op_name in recipes and recipes[op_name].get("kind") == "probe":
        return deep_copy(recipes[op_name])
    raise SuiteError(f"Unknown probe operation: {op_name}")


def normalize_outputs(outputs: Any) -> dict[str, Any]:
    if outputs is None:
        return {}
    if not isinstance(outputs, dict):
        raise SuiteError(f"Outputs must be a mapping, got: {outputs!r}")

    normalized: dict[str, Any] = {}
    for name, spec in outputs.items():
        if isinstance(spec, str):
            normalized[name] = {"discover": {"kind": "path", "path": spec}}
            continue
        if not isinstance(spec, dict):
            raise SuiteError(f"Invalid output specification for {name}: {spec!r}")
        if "discover" in spec:
            normalized[name] = deep_copy(spec)
            continue
        if "path" in spec:
            normalized[name] = {"discover": {"kind": "path", "path": deep_copy(spec["path"])}}
            continue
        if "glob" in spec:
            discover = {"kind": "glob", "pattern": deep_copy(spec["glob"])}
            if "select" in spec:
                discover["select"] = deep_copy(spec["select"])
            normalized[name] = {"discover": discover}
            continue
        raise SuiteError(f"Invalid output specification for {name}: {spec!r}")
    return normalized


def normalize_extract(extract: Any) -> dict[str, Any]:
    if isinstance(extract, str):
        return {"kind": extract}
    if isinstance(extract, dict):
        if "kind" in extract:
            return deep_copy(extract)
        if len(extract) == 1:
            kind, value = next(iter(extract.items()))
            if isinstance(value, dict):
                normalized = {"kind": kind}
                normalized.update(deep_copy(value))
                return normalized
            if kind == "json_path":
                return {"kind": "json_path", "path": deep_copy(value)}
            if kind == "line_pattern":
                return {"kind": "line_pattern", "pattern": deep_copy(value)}
            return {"kind": kind}
    raise SuiteError(f"Invalid extract specification: {extract!r}")


def normalize_recipe(recipe: dict[str, Any]) -> dict[str, Any]:
    normalized = deep_copy(recipe)

    if "probe" in normalized:
        if "kind" in normalized or "command" in normalized:
            raise SuiteError(f"Recipe may not combine 'probe' with 'kind' or 'command': {recipe!r}")
        normalized["kind"] = "probe"
        normalized["command"] = normalized.pop("probe")

    if "action" in normalized:
        if "kind" in normalized or "command" in normalized:
            raise SuiteError(f"Recipe may not combine 'action' with 'kind' or 'command': {recipe!r}")
        normalized["kind"] = "action"
        normalized["command"] = normalized.pop("action")

    defaults = deep_copy(normalized.pop("defaults", {}))
    params = deep_copy(normalized.get("params", {}))
    merged_params = deep_merge(defaults, params)
    if merged_params:
        normalized["params"] = merged_params
    elif "params" in normalized:
        normalized["params"] = {}

    if "cwd" not in normalized and normalized.get("kind") in {"action", "probe"}:
        normalized["cwd"] = "{cwd}"

    if "extract" in normalized:
        normalized["extract"] = normalize_extract(normalized["extract"])

    if "outputs" in normalized:
        normalized["outputs"] = normalize_outputs(normalized["outputs"])

    return normalized


def parse_use_mapping(reference: str) -> dict[str, str]:
    if "." in reference:
        fixture_name, output_name = reference.split(".", 1)
    else:
        fixture_name, output_name = reference, "state"
    if not fixture_name:
        raise SuiteError(f"Invalid fixture output reference: {reference!r}")
    if not output_name:
        raise SuiteError(f"Invalid fixture output reference: {reference!r}")
    return {"fixture": fixture_name, "output": output_name}


def normalize_use(use: Any) -> dict[str, Any]:
    if use is None:
        return {"fixtures": [], "bindings": {}}
    if isinstance(use, str):
        return {"fixtures": [use], "bindings": {}}
    if isinstance(use, list):
        if not all(isinstance(item, str) for item in use):
            raise SuiteError(f"Fixture list must contain only strings: {use!r}")
        return {"fixtures": deep_copy(use), "bindings": {}}
    if isinstance(use, dict):
        bindings: dict[str, dict[str, str]] = {}
        fixtures: list[str] = []
        for alias, reference in use.items():
            if not isinstance(reference, str):
                raise SuiteError(f"Use binding {alias!r} must be a string reference, got: {reference!r}")
            binding = parse_use_mapping(reference)
            bindings[alias] = binding
            if binding["fixture"] not in fixtures:
                fixtures.append(binding["fixture"])
        return {"fixtures": fixtures, "bindings": bindings}
    raise SuiteError(f"Invalid use specification: {use!r}")


def normalize_raw_action(step: dict[str, Any]) -> dict[str, Any]:
    command = deep_copy(step["run"])
    action: dict[str, Any] = {"command": command}
    if "cwd" in step:
        action["cwd"] = deep_copy(step["cwd"])
    outputs = normalize_outputs(step.get("outputs"))
    if outputs:
        action["outputs"] = outputs
    normalized = {"action": action}
    if "save_as" in step:
        normalized["save_as"] = deep_copy(step["save_as"])
    elif outputs:
        normalized["save_as"] = {name: name for name in outputs.keys()}
    return normalized


def normalize_step(step: dict[str, Any], recipes: dict[str, Any]) -> dict[str, Any]:
    action_names = available_action_names(recipes)

    def maybe_rewrite_action(op_name: str, overrides: Any) -> tuple[str, Any]:
        if op_name == "dmrg" and isinstance(overrides, dict) and "sweeps" in overrides:
            return "dmrg_legacy_sweeps", overrides
        return op_name, overrides

    if "run" in step:
        run_spec = step["run"]
        if isinstance(run_spec, (str, list)):
            normalized = normalize_raw_action(step)
        elif isinstance(run_spec, dict):
            run_spec_dict = deep_copy(run_spec)
            if "command" in run_spec_dict:
                pseudo_step = {
                    "run": run_spec_dict["command"],
                    "outputs": run_spec_dict.get("outputs", step.get("outputs")),
                    "save_as": step.get("save_as", run_spec_dict.get("save_as")),
                }
                if "cwd" in run_spec_dict or "cwd" in step:
                    pseudo_step["cwd"] = run_spec_dict.get("cwd", step.get("cwd"))
                normalized = normalize_raw_action(pseudo_step)
            else:
                op_name = run_spec_dict.pop("op", run_spec_dict.pop("recipe", None))
                if op_name is None:
                    raise SuiteError(f"Run step is missing an operation name: {step!r}")
                op_name, overrides = maybe_rewrite_action(op_name, deep_copy(run_spec_dict.get("with", {})))
                normalized = {"action": {"op": op_name, "with": overrides}}
                if "save_as" in step:
                    normalized["save_as"] = deep_copy(step["save_as"])
        else:
            raise SuiteError(f"Invalid run step: {step!r}")
    elif "op" in step or "recipe" in step:
        op_name = step.get("op", step.get("recipe"))
        op_name, overrides = maybe_rewrite_action(op_name, deep_copy(step.get("with", {})))
        normalized = {"action": {"op": op_name, "with": overrides}}
        if "save_as" in step:
            normalized["save_as"] = deep_copy(step["save_as"])
    elif len(step) == 1:
        op_name, overrides = next(iter(step.items()))
        op_name, overrides = maybe_rewrite_action(op_name, overrides)
        if op_name not in action_names:
            raise SuiteError(f"Unknown action shorthand: {op_name}")
        normalized = {
            "action": {
                "op": op_name,
                "with": deep_copy(overrides or {}),
            }
        }
    else:
        raise SuiteError(f"Cannot normalize step: {step!r}")

    if "save_as" not in normalized:
        action_spec = normalized["action"]
        if "op" in action_spec:
            spec = resolve_action_spec(action_spec["op"], recipes)
            outputs = spec.get("outputs", {})
        else:
            outputs = action_spec.get("outputs", {})
        normalized["save_as"] = {name: name for name in outputs.keys()}
    return normalized


def normalize_assertion(assertion: dict[str, Any], recipes: dict[str, Any]) -> dict[str, Any]:
    if "probe" in assertion and "compare" in assertion:
        normalized = deep_copy(assertion)
        probe_spec = normalized["probe"]
        if isinstance(probe_spec, str):
            normalized["probe"] = {"op": probe_spec, "with": deep_copy(assertion.get("with", {}))}
        elif isinstance(probe_spec, dict):
            op_name = probe_spec.pop("op", probe_spec.pop("recipe", None))
            if op_name is None:
                raise SuiteError(f"Probe assertion is missing an operation name: {assertion!r}")
            normalized["probe"] = {"op": op_name, "with": deep_copy(probe_spec.get("with", {}))}
        else:
            raise SuiteError(f"Invalid probe specification: {probe_spec!r}")
        return normalized

    if "probe" in assertion:
        probe_spec = deep_copy(assertion["probe"])
        if isinstance(probe_spec, str):
            probe = {"op": probe_spec, "with": deep_copy(assertion.get("with", {}))}
        elif isinstance(probe_spec, dict):
            op_name = probe_spec.pop("op", probe_spec.pop("recipe", None))
            if op_name is None:
                raise SuiteError(f"Invalid probe specification: {probe_spec!r}")
            probe = {"op": op_name, "with": deep_copy(probe_spec.get("with", {}))}
        else:
            raise SuiteError(f"Invalid probe specification: {probe_spec!r}")
        compare = build_compare(assertion)
        return {"probe": probe, "compare": compare}

    if len(assertion) == 1:
        op_name, spec = next(iter(assertion.items()))
        if op_name in available_probe_names(recipes):
            spec_dict = deep_copy(spec or {})
            compare = build_compare(spec_dict)
            probe = {"op": op_name, "with": strip_compare_fields(spec_dict)}
            return {"probe": probe, "compare": compare}

    raise SuiteError(f"Cannot normalize assertion: {assertion!r}")


def normalize_fixture(fixture: dict[str, Any], recipes: dict[str, Any]) -> dict[str, Any]:
    normalized = deep_copy(fixture)
    if "run" in normalized and "steps" not in normalized and "build" not in normalized:
        build_spec: dict[str, Any] = {"run": normalized.pop("run")}
        if "outputs" in normalized:
            build_spec["outputs"] = normalized.pop("outputs")
        if "save_as" in normalized:
            build_spec["save_as"] = normalized.pop("save_as")
        if "cwd" in normalized:
            build_spec["cwd"] = normalized.pop("cwd")
        normalized["steps"] = [build_spec]
    if "build" in normalized and "steps" not in normalized:
        build_spec = normalized.pop("build")
        if isinstance(build_spec, str):
            normalized["steps"] = [{"recipe": build_spec}]
        elif isinstance(build_spec, dict):
            normalized["steps"] = [build_spec]
        else:
            raise SuiteError(f"Invalid fixture build specification: {build_spec!r}")
    if "steps" in normalized:
        normalized["steps"] = [normalize_step(step, recipes) for step in as_list(normalized["steps"])]
    if "certify" in normalized:
        normalized["certify"] = [
            normalize_assertion(assertion, recipes)
            for assertion in as_list(normalized["certify"])
        ]
    return normalized


def normalize_test(test: dict[str, Any], recipes: dict[str, Any]) -> dict[str, Any]:
    normalized = deep_copy(test)
    if "use" in normalized:
        normalized["use"] = normalize_use(normalized["use"])
    else:
        normalized["use"] = normalize_use(None)
    if "steps" in normalized:
        normalized["steps"] = [normalize_step(step, recipes) for step in as_list(normalized["steps"])]
    if "expect" in normalized:
        normalized["expect"] = [
            normalize_assertion(assertion, recipes)
            for assertion in as_list(normalized["expect"])
        ]
    return normalized


def normalize_suite(suite: dict[str, Any]) -> dict[str, Any]:
    normalized = deep_merge(DEFAULT_OVERLAY, deep_copy(suite))
    recipes = {
        name: normalize_recipe(spec)
        for name, spec in normalized.get("recipes", {}).items()
    }
    normalized["recipes"] = recipes
    normalized["fixtures"] = {
        name: normalize_fixture(spec, recipes)
        for name, spec in normalized.get("fixtures", {}).items()
    }
    normalized["tests"] = {
        name: normalize_test(spec, recipes)
        for name, spec in normalized.get("tests", {}).items()
    }
    return normalized


class SuiteRunner:
    def __init__(
        self,
        suite: dict[str, Any],
        bin_dir: Path,
        work_root: Path,
        config: str | None = None,
        verbose: bool = False,
        explain: bool = False,
        trace: bool = False,
    ):
        self.suite = suite
        self.bin_dir = bin_dir.resolve()
        self.config = normalize_config(config)
        self.bin_path_dirs = executable_search_dirs(self.bin_dir, self.config)
        self.repo_root = Path(__file__).resolve().parents[1]
        self.work_root = work_root
        self.verbose = verbose
        self.explain = explain
        self.trace = trace
        self.defaults = suite.get("defaults", {})
        self.recipes = suite.get("recipes", {})
        self.fixtures = suite.get("fixtures", {})
        self.tests = suite.get("tests", {})
        self.fixture_root = self.work_root / "fixture-cache"
        self.run_root = self.work_root / "runs"
        ensure_directory(self.fixture_root)
        ensure_directory(self.run_root)
        self.built_fixtures: dict[str, dict[str, Any]] = {}
        self.validate_fixtures()

    def log(self, message: str) -> None:
        if self.verbose:
            print(message)

    def log_explain(self, message: str) -> None:
        if self.explain:
            print(message)

    def log_trace(self, message: str) -> None:
        if self.trace:
            print(message)

    def make_base_context(self) -> dict[str, Any]:
        return {
            "bin_dir": str(self.bin_dir),
            "config": self.config or "",
            "repo_root": str(self.repo_root),
            "work_root": str(self.work_root),
        }

    def merged_env(self, recipe: dict[str, Any], context: dict[str, Any]) -> dict[str, str]:
        env = dict(os.environ)
        default_env = render_value(deep_copy(self.defaults.get("env", {})), context)
        recipe_env = render_value(deep_copy(recipe.get("env", {})), context)
        env.update({key: str(value) for key, value in default_env.items()})
        env.update({key: str(value) for key, value in recipe_env.items()})
        bin_path_entries = [
            str(path.resolve())
            for path in self.bin_path_dirs
            if path.is_dir()
        ]
        if bin_path_entries:
            existing_path = env.get("PATH", "")
            env["PATH"] = os.pathsep.join(bin_path_entries + [existing_path])
        if "MP_BINPATH" in env:
            ensure_directory(Path(env["MP_BINPATH"]))
        return env

    def resolve_cwd(self, recipe: dict[str, Any], context: dict[str, Any]) -> Path:
        cwd_value = recipe.get("cwd", "{scratch}")
        cwd_text = render_value(cwd_value, context)
        cwd_path = Path(cwd_text)
        if not cwd_path.is_absolute():
            cwd_path = Path(context["scratch"]) / cwd_path
        ensure_directory(cwd_path)
        return cwd_path

    def render_command(self, recipe: dict[str, Any], call_context: dict[str, Any]) -> list[str]:
        command = deep_copy(recipe["command"])
        if isinstance(command, str):
            argv = [render_value(token, call_context) for token in shlex.split(command)]
        elif isinstance(command, list) and all(isinstance(item, str) for item in command):
            argv = [render_value(item, call_context) for item in command]
        else:
            raise SuiteError("Command must render to a string or a list of strings")
        if argv and "/" not in argv[0]:
            candidate = self.resolve_executable(argv[0])
            if candidate is not None:
                argv = [str(candidate)] + argv[1:]
        return argv

    def resolve_executable(self, program: str):
        for directory in self.bin_path_dirs:
            candidate = directory / program
            if candidate.exists() and os.access(candidate, os.X_OK):
                return candidate.resolve()
        return None

    def render_command_text(self, command: list[str]) -> str:
        return shlex.join(command)

    def resolve_action_spec(self, op_name: str) -> dict[str, Any]:
        return resolve_action_spec(op_name, self.recipes)

    def resolve_probe_spec(self, op_name: str) -> dict[str, Any]:
        return resolve_probe_spec(op_name, self.recipes)

    def validate_fixtures(self) -> None:
        visiting: list[str] = []
        visited: set[str] = set()

        def dfs(name: str) -> None:
            if name in visited:
                return
            if name in visiting:
                cycle = visiting[visiting.index(name):] + [name]
                raise SuiteError(f"Fixture dependency cycle detected: {' -> '.join(cycle)}")
            if name not in self.fixtures:
                raise SuiteError(f"Unknown fixture dependency: {name}")

            visiting.append(name)
            parent = self.fixtures[name].get("from")
            if parent is not None:
                if parent not in self.fixtures:
                    raise SuiteError(f"Fixture {name} depends on unknown fixture {parent}")
                dfs(parent)
            visiting.pop()
            visited.add(name)

        for fixture_name in self.fixtures:
            dfs(fixture_name)

    def materialize_action(self, action_spec: dict[str, Any]) -> tuple[str, dict[str, Any], dict[str, Any]]:
        if "command" in action_spec:
            command = action_spec["command"]
            if isinstance(command, str):
                action_name = shlex.split(command)[0] if command.strip() else "command"
            elif isinstance(command, list) and command:
                action_name = str(command[0])
            else:
                action_name = "command"
            runtime_spec = {
                "command": deep_copy(command),
                "cwd": deep_copy(action_spec.get("cwd", "{cwd}")),
                "outputs": deep_copy(action_spec.get("outputs", {})),
            }
            return (action_name, runtime_spec, {})
        if "op" in action_spec:
            return (
                action_spec["op"],
                self.resolve_action_spec(action_spec["op"]),
                deep_copy(action_spec.get("with", {})),
            )
        raise SuiteError(f"Invalid action specification: {action_spec!r}")

    def unexpected_created_paths(
        self,
        cwd: Path,
        before: set[Path],
        after: set[Path],
        discovered: dict[str, str],
    ) -> list[Path]:
        cwd = cwd.resolve()
        new_paths = sorted(after - before)
        if not new_paths:
            return []

        declared = [Path(path_text).resolve() for path_text in discovered.values()]

        def is_expected(relative_path: Path) -> bool:
            absolute_path = (cwd / relative_path).resolve()
            for output_path in declared:
                if absolute_path == output_path:
                    return True
                if is_relative_to(output_path, absolute_path):
                    return True
                if output_path.exists() and output_path.is_dir() and is_relative_to(absolute_path, output_path):
                    return True
            return False

        return [path for path in new_paths if not is_expected(path)]

    def apply_scope_defaults(self, overrides: dict[str, Any], context: dict[str, Any]) -> dict[str, Any]:
        scoped = deep_copy(overrides)
        fixture_name = scoped.pop("fixture", None)
        if fixture_name is not None:
            fixtures = context.get("fixtures", {})
            if fixture_name not in fixtures:
                raise SuiteError(f"Unknown fixture scope: {fixture_name}")
            default_cwd = fixtures[fixture_name]["dir"]
        else:
            default_cwd = context.get("_default_cwd", context["scratch"])

        if "cwd" not in scoped and default_cwd is not None:
            scoped["cwd"] = default_cwd
        bound_paths = context.get("bindings", {})
        current_state_path = context.get("_current_state_path")
        for key in ("state", "psi1", "psi2"):
            if key in scoped:
                continue
            if key in bound_paths:
                scoped[key] = bound_paths[key]
                continue
            if current_state_path is not None:
                scoped[key] = current_state_path
        return scoped

    def discover_outputs(
        self,
        recipe: dict[str, Any],
        call_context: dict[str, Any],
        cwd: Path,
    ) -> dict[str, str]:
        discovered: dict[str, str] = {}
        for output_name, output_spec in recipe.get("outputs", {}).items():
            discover = output_spec.get("discover", {})
            kind = discover.get("kind")
            if kind == "path":
                rendered_path = render_value(discover["path"], call_context)
                path = Path(rendered_path)
                if not path.is_absolute():
                    path = cwd / path
                if not path.exists():
                    raise SuiteError(f"Declared output path was not created: {path}")
                discovered[output_name] = str(path.resolve())
                continue
            if kind == "glob":
                rendered_pattern = render_value(discover["pattern"], call_context)
                pattern_path = Path(rendered_pattern)
                if not pattern_path.is_absolute():
                    pattern_path = cwd / pattern_path
                matches = sorted(pattern_path.parent.glob(pattern_path.name))
                if not matches:
                    raise SuiteError(f"No outputs matched glob {pattern_path}")
                select = discover.get("select", "first")
                if select == "newest":
                    chosen = max(matches, key=lambda path: (path.stat().st_mtime_ns, str(path)))
                else:
                    chosen = matches[0]
                discovered[output_name] = str(chosen.resolve())
                continue
            raise SuiteError(f"Unsupported output discovery kind: {kind}")
        return discovered

    def describe_outputs(
        self,
        recipe: dict[str, Any],
        call_context: dict[str, Any],
        cwd: Path,
    ) -> list[str]:
        descriptions: list[str] = []
        for output_name, output_spec in recipe.get("outputs", {}).items():
            discover = output_spec.get("discover", {})
            kind = discover.get("kind")
            if kind == "path":
                rendered_path = render_value(discover["path"], call_context)
                path = Path(rendered_path)
                if not path.is_absolute():
                    path = cwd / path
                descriptions.append(f"{output_name}: {path.resolve()}")
                continue
            if kind == "glob":
                rendered_pattern = render_value(discover["pattern"], call_context)
                pattern_path = Path(rendered_pattern)
                if not pattern_path.is_absolute():
                    pattern_path = cwd / pattern_path
                descriptions.append(
                    f"{output_name}: {pattern_path.resolve()} (select={discover.get('select', 'first')})"
                )
                continue
            descriptions.append(f"{output_name}: <unsupported discovery {kind}>")
        return descriptions

    def run_action(self, action_spec: dict[str, Any], context: dict[str, Any]) -> dict[str, str]:
        action_name, spec, overrides = self.materialize_action(action_spec)
        params = deep_copy(spec.get("params", {}))
        params.update(self.apply_scope_defaults(overrides, context))
        call_context = deep_copy(context)
        call_context.update(render_value(params, context))
        cwd = self.resolve_cwd(spec, call_context)
        env = self.merged_env(spec, call_context)
        command = self.render_command(spec, call_context)
        ignored_roots: list[Path] = []
        if "MP_BINPATH" in env:
            mp_binpath = Path(env["MP_BINPATH"]).resolve()
            if is_relative_to(mp_binpath, cwd.resolve()):
                ignored_roots.append(mp_binpath)
        before_snapshot = snapshot_tree(cwd, ignore_roots=ignored_roots)
        self.log(f"[action] {' '.join(command)}")
        if self.explain:
            print(f"ACTION {action_name}")
            print(f"  cwd: {cwd}")
            print(f"  command: {self.render_command_text(command)}")
            outputs = self.describe_outputs(spec, call_context, cwd)
            if outputs:
                print("  outputs:")
                for output in outputs:
                    print(f"    {output}")
        result = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if self.trace:
            print(f"TRACE action {action_name} exit={result.returncode}")
            print("  stdout:")
            print(indent_block(result.stdout, "    "))
            print("  stderr:")
            print(indent_block(result.stderr, "    "))
        if result.returncode != 0:
            raise SuiteError(
                f"Action {action_name} failed with exit code {result.returncode}\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )
        discovered = self.discover_outputs(spec, call_context, cwd)
        after_snapshot = snapshot_tree(cwd, ignore_roots=ignored_roots)
        unexpected = self.unexpected_created_paths(cwd, before_snapshot, after_snapshot, discovered)
        if unexpected:
            unexpected_text = "\n".join(f"  {path}" for path in unexpected)
            raise SuiteError(
                f"Action {action_name} created undeclared outputs in {cwd}:\n{unexpected_text}"
            )
        if self.trace and discovered:
            print("  discovered outputs:")
            for name, value in discovered.items():
                print(f"    {name}: {value}")
        return discovered

    def run_probe(
        self,
        op_name: str,
        overrides: dict[str, Any],
        context: dict[str, Any],
    ) -> Any:
        spec = self.resolve_probe_spec(op_name)
        params = deep_copy(spec.get("params", {}))
        params.update(self.apply_scope_defaults(overrides, context))
        call_context = deep_copy(context)
        call_context.update(render_value(params, context))
        cwd = self.resolve_cwd(spec, call_context)
        env = self.merged_env(spec, call_context)
        command = self.render_command(spec, call_context)
        self.log(f"[probe] {' '.join(command)}")
        if self.explain:
            print(f"PROBE {op_name}")
            print(f"  cwd: {cwd}")
            print(f"  command: {self.render_command_text(command)}")
        result = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if self.trace:
            print(f"TRACE probe {op_name} exit={result.returncode}")
            print("  stdout:")
            print(indent_block(result.stdout, "    "))
            print("  stderr:")
            print(indent_block(result.stderr, "    "))
        if result.returncode != 0:
            raise SuiteError(
                f"Probe {op_name} failed with exit code {result.returncode}\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )
        extract = render_value(deep_copy(spec.get("extract", {})), call_context)
        kind = extract.get("kind")
        if kind == "last_float":
            value = parse_float(last_nonblank_line(result.stdout))
        elif kind == "last_complex":
            value = parse_complex(last_nonblank_line(result.stdout))
        elif kind == "last_text":
            value = last_nonblank_line(result.stdout)
        elif kind == "full_text":
            value = normalize_capture(result.stdout)
        elif kind == "json_path":
            payload = json.loads(normalize_capture(result.stdout))
            value = extract_json_path(payload, extract.get("path", ""))
        elif kind == "line_pattern":
            value = extract_line_pattern(
                result.stdout,
                str(extract.get("pattern", "")),
                extract.get("capture"),
            )
        else:
            raise SuiteError(f"Unsupported probe extractor kind: {kind}")
        if self.trace:
            print(f"  extracted value: {value!r}")
        return value

    def evaluate_assertion(self, assertion: dict[str, Any], context: dict[str, Any]) -> None:
        probe_spec = assertion["probe"]
        if self.explain:
            print(f"ASSERT {probe_spec['op']}")
        actual = self.run_probe(probe_spec["op"], probe_spec.get("with", {}), context)
        compare = assertion["compare"]
        kind = compare["kind"]
        if kind == "approx":
            expected = coerce_expected(compare["value"], actual)
            abs_tol = float(compare.get("abs", 0.0))
            rel_tol = float(compare.get("rel", 0.0))
            if not approx_equal(actual, expected, abs_tol=abs_tol, rel_tol=rel_tol):
                raise SuiteError(
                    f"Approx comparison failed: actual={actual!r} expected={expected!r} "
                    f"abs={abs_tol} rel={rel_tol}"
                )
            if self.trace:
                print(
                    f"  assertion passed: approx actual={actual!r} expected={expected!r} "
                    f"abs={abs_tol} rel={rel_tol}"
                )
            return
        if kind == "equals":
            expected = coerce_expected(compare["value"], actual)
            if actual != expected:
                raise SuiteError(f"Equality comparison failed: actual={actual!r} expected={expected!r}")
            if self.trace:
                print(f"  assertion passed: equals actual={actual!r} expected={expected!r}")
            return
        if kind == "matches":
            pattern = str(compare["value"])
            if re.search(pattern, str(actual)) is None:
                raise SuiteError(f"Regex comparison failed: actual={actual!r} pattern={pattern!r}")
            if self.trace:
                print(f"  assertion passed: matches actual={actual!r} pattern={pattern!r}")
            return
        if kind == "contains":
            expected = compare["value"]
            if isinstance(actual, str):
                ok = str(expected) in actual
            else:
                try:
                    ok = expected in actual
                except TypeError:
                    ok = False
            if not ok:
                raise SuiteError(f"Contains comparison failed: actual={actual!r} expected={expected!r}")
            if self.trace:
                print(f"  assertion passed: contains actual={actual!r} expected={expected!r}")
            return
        if kind == "exists":
            should_exist = bool(compare.get("value", True))
            if isinstance(actual, (str, os.PathLike)):
                exists = Path(actual).exists()
            else:
                exists = actual is not None
            if exists != should_exist:
                raise SuiteError(
                    f"Exists comparison failed: actual={actual!r} exists={exists!r} expected={should_exist!r}"
                )
            if self.trace:
                print(
                    f"  assertion passed: exists actual={actual!r} exists={exists!r} expected={should_exist!r}"
                )
            return
        raise SuiteError(f"Unsupported comparison kind: {kind}")

    def run_steps(
        self,
        step_list: list[dict[str, Any]],
        context: dict[str, Any],
        output_namespace: dict[str, dict[str, str]],
    ) -> dict[str, str]:
        published: dict[str, str] = {}
        for step in step_list:
            action_spec = step["action"]
            outputs = self.run_action(action_spec, context)
            for output_name, alias in step.get("save_as", {}).items():
                if output_name not in outputs:
                    raise SuiteError(f"Step tried to save missing output {output_name}")
                output_namespace.setdefault(output_name, {})
                output_namespace[output_name][alias] = outputs[output_name]
                published[alias] = outputs[output_name]
                if output_name == "state":
                    context["_current_state_path"] = outputs[output_name]
                    context["bindings"]["state"] = outputs[output_name]
                context["bindings"][alias] = outputs[output_name]
            context["artifacts"].update(published)
        return published

    def fixture_steps(self, fixture_spec: dict[str, Any]) -> list[dict[str, Any]]:
        if "steps" in fixture_spec:
            return fixture_spec["steps"]
        if "build" in fixture_spec:
            return [{"action": fixture_spec["build"], "save_as": {}}]
        raise SuiteError("Fixture must define steps or build")

    def inherit_fixture_outputs(
        self,
        source_info: dict[str, Any] | None,
        target_dir: Path,
    ) -> dict[str, str]:
        if source_info is None:
            return {}
        inherited: dict[str, str] = {}
        source_dir = Path(source_info["dir"])
        for name, path_text in source_info.items():
            if name == "dir":
                continue
            path = Path(path_text)
            try:
                relative = path.relative_to(source_dir)
            except ValueError:
                inherited[name] = path_text
                continue
            inherited[name] = str((target_dir / relative).resolve())
        return inherited

    def materialize_fixture(self, fixture_name: str, fixture_info: dict[str, Any], scratch: Path) -> dict[str, Any]:
        fixture_root = scratch / "fixtures"
        ensure_directory(fixture_root)
        fixture_dir = fixture_root / fixture_name
        shutil.copytree(fixture_info["dir"], fixture_dir)
        local_info: dict[str, Any] = {"dir": str(fixture_dir.resolve())}
        local_info.update(self.inherit_fixture_outputs(fixture_info, fixture_dir.resolve()))
        return local_info

    def build_fixture(self, name: str) -> dict[str, Any]:
        if name in self.built_fixtures:
            self.log_explain(f"FIXTURE {name} (cached)")
            return self.built_fixtures[name]

        fixture_spec = self.fixtures[name]
        fixture_dir = self.fixture_root / name
        self.log_explain(f"FIXTURE {name}")
        if "from" in fixture_spec:
            self.log_explain(f"  depends on: {fixture_spec['from']}")
        self.log_explain(f"  cache dir: {fixture_dir}")
        if fixture_dir.exists():
            shutil.rmtree(fixture_dir)

        parent_info: dict[str, Any] | None = None
        if "from" in fixture_spec:
            parent_info = self.build_fixture(fixture_spec["from"])
            shutil.copytree(parent_info["dir"], fixture_dir)
        else:
            ensure_directory(fixture_dir)

        inherited_outputs = self.inherit_fixture_outputs(parent_info, fixture_dir)
        parent_context = {"dir": str(fixture_dir), **inherited_outputs} if parent_info is not None else {}
        bindings: dict[str, str] = {}
        if "state" in inherited_outputs:
            bindings["state"] = inherited_outputs["state"]

        base_context = self.make_base_context()
        context: dict[str, Any] = {
            **base_context,
            "scratch": str(fixture_dir),
            "_default_cwd": str(fixture_dir),
            "_current_state_path": bindings.get("state"),
            "bindings": deep_copy(bindings),
            "fixture": {"dir": str(fixture_dir)},
            "parent": deep_copy(parent_context),
            "fixtures": deep_copy(self.built_fixtures),
            "artifacts": deep_copy(bindings),
        }
        output_namespace: dict[str, dict[str, str]] = {}
        published = self.run_steps(self.fixture_steps(fixture_spec), context, output_namespace)
        fixture_info: dict[str, Any] = {"dir": str(fixture_dir)}
        fixture_info.update(inherited_outputs)
        fixture_info.update(published)
        if "state" in fixture_info:
            context["_current_state_path"] = fixture_info["state"]
            context["bindings"]["state"] = fixture_info["state"]
        context["fixture"].update(fixture_info)
        if self.trace and fixture_info:
            print("  fixture outputs:")
            for key, value in fixture_info.items():
                print(f"    {key}: {value}")

        for assertion in fixture_spec.get("certify", []):
            self.evaluate_assertion(assertion, context)

        self.built_fixtures[name] = fixture_info
        return fixture_info

    def run_test(self, name: str) -> None:
        test_spec = self.tests[name]
        use_spec = test_spec.get("use", {"fixtures": [], "bindings": {}})
        self.log_explain(f"TEST {name}")
        if use_spec.get("fixtures"):
            self.log_explain(f"  uses fixtures: {', '.join(use_spec['fixtures'])}")
        if use_spec.get("bindings"):
            self.log_explain("  bindings:")
            for alias, binding in use_spec["bindings"].items():
                self.log_explain(f"    {alias}: {binding['fixture']}.{binding['output']}")
        scratch = self.run_root / name
        if scratch.exists():
            shutil.rmtree(scratch)
        ensure_directory(scratch)
        self.log_explain(f"  scratch dir: {scratch}")

        fixture_cache_infos: dict[str, Any] = {}
        for fixture_name in use_spec.get("fixtures", []):
            fixture_cache_infos[fixture_name] = self.build_fixture(fixture_name)

        fixtures: dict[str, Any] = {}
        for fixture_name, fixture_info in fixture_cache_infos.items():
            fixtures[fixture_name] = self.materialize_fixture(fixture_name, fixture_info, scratch)
            self.log_explain(
                f"  materialize fixture {fixture_name} -> {fixtures[fixture_name]['dir']}"
            )

        base_context = self.make_base_context()
        bindings: dict[str, str] = {}
        for alias, binding in use_spec.get("bindings", {}).items():
            fixture_name = binding["fixture"]
            output_name = binding["output"]
            if fixture_name not in fixtures:
                raise SuiteError(f"Use binding references unknown fixture {fixture_name}")
            if output_name not in fixtures[fixture_name]:
                raise SuiteError(
                    f"Use binding {alias} references missing output {fixture_name}.{output_name}"
                )
            bindings[alias] = fixtures[fixture_name][output_name]

        if not bindings and len(fixtures) == 1:
            only_fixture = next(iter(fixtures.values()))
            if "state" in only_fixture:
                bindings["state"] = only_fixture["state"]

        default_cwd = str(scratch)
        if len(fixtures) == 1:
            default_cwd = next(iter(fixtures.values()))["dir"]
        elif use_spec.get("bindings"):
            bound_fixture_names = {binding["fixture"] for binding in use_spec["bindings"].values()}
            if len(bound_fixture_names) == 1:
                default_cwd = fixtures[next(iter(bound_fixture_names))]["dir"]

        current_state_path = bindings.get("state")
        context: dict[str, Any] = {
            **base_context,
            "scratch": str(scratch),
            "_default_cwd": default_cwd,
            "_current_state_path": current_state_path,
            "fixtures": deep_copy(fixtures),
            "bindings": deep_copy(bindings),
            "state": {},
            "lattice": {},
            "artifacts": deep_copy(bindings),
        }

        output_namespace: dict[str, dict[str, str]] = {"state": context["state"], "lattice": context["lattice"]}
        self.run_steps(test_spec.get("steps", []), context, output_namespace)

        for assertion in test_spec.get("expect", []):
            self.evaluate_assertion(assertion, context)

    def run(self, selected_tests: list[str] | None = None) -> int:
        test_names = selected_tests or list(self.tests.keys())
        unknown = [name for name in test_names if name not in self.tests]
        if unknown:
            raise SuiteError(f"Unknown test(s): {', '.join(unknown)}")
        failures: list[tuple[str, str]] = []
        for name in test_names:
            try:
                self.run_test(name)
                print(f"PASS {name}")
            except Exception as exc:  # pragma: no cover - suite failure path
                failures.append((name, str(exc)))
                print(f"FAIL {name}")
                print(str(exc))
        passed = len(test_names) - len(failures)
        print(f"Summary: {len(test_names)} total, {passed} passed, {len(failures)} failed")
        if failures:
            print("Failed tests:")
            for name, _message in failures:
                print(f"- {name}")
            return 1
        return 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MPToolkit integration test suite")
    parser.add_argument("suite", type=Path, help="Path to YAML suite file")
    parser.add_argument("--bin-dir", type=Path, required=True, help="Directory containing MPToolkit binaries")
    parser.add_argument("--config", help="Prefer binaries from this CMake configuration when searching a multi-config build tree")
    parser.add_argument("--work-root", type=Path, help="Working directory for fixture and test outputs")
    parser.add_argument("--verbose", action="store_true", help="Print executed commands")
    parser.add_argument("--explain", action="store_true", help="Explain fixture dependencies, files, and resolved commands")
    parser.add_argument("--trace", action="store_true", help="Trace command output, extracted values, and assertion results")
    parser.add_argument("--dump-ir", action="store_true", help="Dump the normalized internal suite representation and exit")
    parser.add_argument("tests", nargs="*", help="Optional list of specific tests to run")
    if hasattr(parser, "parse_intermixed_args"):
        return parser.parse_intermixed_args(argv)
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    try:
        args = parse_args(argv)
        suite = normalize_suite(load_raw_suite(args.suite.resolve()))
        if args.dump_ir:
            yaml.safe_dump(suite, sys.stdout, sort_keys=False)
            return 0

        def run_with_work_root(work_root: Path) -> int:
            ensure_directory(work_root)

            runner = SuiteRunner(
                suite=suite,
                bin_dir=args.bin_dir,
                work_root=work_root,
                config=args.config,
                verbose=args.verbose,
                explain=args.explain,
                trace=args.trace,
            )
            return runner.run(selected_tests=args.tests or None)

        if args.work_root is None:
            with tempfile.TemporaryDirectory(prefix="mptk-tests-") as temp_root:
                return run_with_work_root(Path(temp_root))

        return run_with_work_root(args.work_root)
    except SuiteError as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
