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

import yaml


ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
TEMPLATE_RE = re.compile(r"\{([^{}]+)\}")
DEFAULT_OVERLAY = {
    "defaults": {
        "env": {
            "LC_ALL": "C",
            "OMP_NUM_THREADS": "1",
            "MP_BINPATH": "{scratch}/tmp",
        }
    }
}


class SuiteError(RuntimeError):
    pass


def strip_ansi(text: str) -> str:
    return ANSI_ESCAPE_RE.sub("", text)


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def deep_copy(value: Any) -> Any:
    return copy.deepcopy(value)


def resolve_reference(expr: str, context: dict[str, Any]) -> Any:
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
        if isinstance(probe_spec, dict):
            op_name = probe_spec.pop("op", probe_spec.pop("recipe", None))
            if op_name is None:
                raise SuiteError(f"Probe assertion is missing an operation name: {assertion!r}")
            normalized["probe"] = {"op": op_name, "with": deep_copy(probe_spec.get("with", {}))}
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
        normalized["use"] = as_list(normalized["use"])
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
        verbose: bool = False,
        explain: bool = False,
        trace: bool = False,
    ):
        self.suite = suite
        self.bin_dir = bin_dir
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
            "work_root": str(self.work_root),
        }

    def merged_env(self, recipe: dict[str, Any], context: dict[str, Any]) -> dict[str, str]:
        env = dict(os.environ)
        default_env = render_value(deep_copy(self.defaults.get("env", {})), context)
        recipe_env = render_value(deep_copy(recipe.get("env", {})), context)
        env.update({key: str(value) for key, value in default_env.items()})
        env.update({key: str(value) for key, value in recipe_env.items()})
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
        command = render_value(deep_copy(recipe["command"]), call_context)
        if isinstance(command, str):
            argv = shlex.split(command)
        elif isinstance(command, list) and all(isinstance(item, str) for item in command):
            argv = command
        else:
            raise SuiteError("Command must render to a string or a list of strings")
        if argv and "/" not in argv[0]:
            candidate = self.bin_dir / argv[0]
            if candidate.exists():
                argv = [str(candidate)] + argv[1:]
        return argv

    def render_command_text(self, command: list[str]) -> str:
        return shlex.join(command)

    def resolve_action_spec(self, op_name: str) -> dict[str, Any]:
        return resolve_action_spec(op_name, self.recipes)

    def resolve_probe_spec(self, op_name: str) -> dict[str, Any]:
        return resolve_probe_spec(op_name, self.recipes)

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
        default_state_name = context.get("_default_state_name")
        if default_state_name is not None:
            for key in ("state", "psi1", "psi2"):
                if key not in scoped:
                    scoped[key] = default_state_name
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
                    chosen = max(matches, key=lambda path: path.stat().st_mtime)
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
        elif kind == "json_path":
            payload = json.loads(normalize_capture(result.stdout))
            value = extract_json_path(payload, extract.get("path", ""))
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
                    context["_default_state_name"] = Path(outputs[output_name]).name
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

        base_context = self.make_base_context()
        context: dict[str, Any] = {
            **base_context,
            "scratch": str(fixture_dir),
            "_default_cwd": str(fixture_dir),
            "_default_state_name": None,
            "fixture": {"dir": str(fixture_dir)},
            "parent": deep_copy(parent_context),
            "fixtures": deep_copy(self.built_fixtures),
            "artifacts": {},
        }
        output_namespace: dict[str, dict[str, str]] = {}
        published = self.run_steps(self.fixture_steps(fixture_spec), context, output_namespace)
        fixture_info: dict[str, Any] = {"dir": str(fixture_dir)}
        fixture_info.update(inherited_outputs)
        fixture_info.update(published)
        if "state" in fixture_info:
            context["_default_state_name"] = Path(fixture_info["state"]).name
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
        self.log_explain(f"TEST {name}")
        if test_spec.get("use"):
            self.log_explain(f"  uses fixtures: {', '.join(test_spec['use'])}")
        scratch = self.run_root / name
        if scratch.exists():
            shutil.rmtree(scratch)
        ensure_directory(scratch)
        self.log_explain(f"  scratch dir: {scratch}")

        fixture_cache_infos: dict[str, Any] = {}
        for fixture_name in test_spec.get("use", []):
            fixture_cache_infos[fixture_name] = self.build_fixture(fixture_name)

        fixtures: dict[str, Any] = {}
        for fixture_name, fixture_info in fixture_cache_infos.items():
            fixtures[fixture_name] = self.materialize_fixture(fixture_name, fixture_info, scratch)
            self.log_explain(
                f"  materialize fixture {fixture_name} -> {fixtures[fixture_name]['dir']}"
            )

        base_context = self.make_base_context()
        default_cwd = None
        default_state_name = None
        if len(fixtures) == 1:
            only_fixture = next(iter(fixtures.values()))
            default_cwd = only_fixture["dir"]
            if "state" in only_fixture:
                default_state_name = Path(only_fixture["state"]).name
        context: dict[str, Any] = {
            **base_context,
            "scratch": str(scratch),
            "_default_cwd": default_cwd,
            "_default_state_name": default_state_name,
            "fixtures": deep_copy(fixtures),
            "state": {},
            "lattice": {},
            "artifacts": {},
        }

        output_namespace: dict[str, dict[str, str]] = {"state": context["state"], "lattice": context["lattice"]}
        self.run_steps(test_spec.get("steps", []), context, output_namespace)

        for assertion in test_spec.get("expect", []):
            self.evaluate_assertion(assertion, context)

    def run(self, selected_tests: list[str] | None = None) -> int:
        test_names = selected_tests or list(self.tests.keys())
        failures: list[tuple[str, str]] = []
        for name in test_names:
            try:
                self.run_test(name)
                print(f"PASS {name}")
            except Exception as exc:  # pragma: no cover - prototype failure path
                failures.append((name, str(exc)))
                print(f"FAIL {name}")
                print(str(exc))
        if failures:
            print(f"{len(failures)} test(s) failed")
            return 1
        print(f"{len(test_names)} test(s) passed")
        return 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MPToolkit declarative test prototype")
    parser.add_argument("suite", type=Path, help="Path to YAML suite file")
    parser.add_argument("--bin-dir", type=Path, required=True, help="Directory containing MPToolkit binaries")
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
    args = parse_args(argv)
    suite = normalize_suite(load_raw_suite(args.suite.resolve()))
    if args.dump_ir:
        yaml.safe_dump(suite, sys.stdout, sort_keys=False)
        return 0

    work_root = args.work_root
    if work_root is None:
        work_root = Path(tempfile.mkdtemp(prefix="mptk-test-prototype-"))
    ensure_directory(work_root)

    runner = SuiteRunner(
        suite=suite,
        bin_dir=args.bin_dir,
        work_root=work_root,
        verbose=args.verbose,
        explain=args.explain,
        trace=args.trace,
    )
    return runner.run(selected_tests=args.tests or None)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
