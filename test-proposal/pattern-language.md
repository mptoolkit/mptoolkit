# Pattern Language

## Purpose

This document defines the default text-matching semantics for probe recipes and
text-based expectations.

The main goal is to let test authors match useful output without writing raw
regex for common cases.

## Default Normalization

Unless explicitly overridden, text matching should:

- read UTF-8 text from the selected stream
- strip ANSI escape sequences
- normalize line endings
- trim leading and trailing whitespace per line
- use flexible horizontal whitespace matching

Flexible whitespace means that a literal space in the pattern matches zero or
more spaces or tabs in the output.

The current prototype uses `%(...)` for typed captures so that ordinary suite
template interpolation can continue to use `{...}` without ambiguity.

Example:

```yaml
stdout: "Energy = %(float)"
```

should match all of:

```text
Energy = -0.443
Energy=-0.443
Energy    =    -0.443
```

## Matching Modes

Recommended built-in selectors:

- `stdout`: find a matching line anywhere in stdout
- `stderr`: find a matching line anywhere in stderr
- `merged`: find a matching line in merged stdout/stderr
- `stdout_starts_with`
- `stdout_contains`
- `stdout_last`
- `stdout_first`

Examples:

```yaml
stdout: "Energy = %(float)"
```

```yaml
stdout_starts_with: "Converged after %(int) sweeps"
```

```yaml
stdout_contains: "Lanczos terminated early"
```

```yaml
stdout_last: "%(float)"
```

## Typed Placeholders

Recommended standard placeholders:

- `%(float)`
- `%(int)`
- `%(complex)`
- `%(word)`
- `%(rest)`

Named captures should also be supported:

- `%(float:energy)`
- `%(int:steps)`
- `%(word:outfile)`

The engine should compile these placeholders to robust internal regex and then
parse them into typed values.

Test authors should not need to write numeric regex directly for ordinary use.

## Pattern Semantics

Literal text matches literally after whitespace normalization.

Punctuation should be matched permissively with optional surrounding whitespace
where that helps common CLI output.

For example:

```yaml
stdout: "Energy = %(float)"
```

should match `Energy=-0.2` as well as `Energy = -0.2`.

Patterns are line-oriented by default. A pattern matches one line unless a
recipe explicitly opts into multi-line capture.

## Extraction

If a pattern contains exactly one unnamed typed placeholder, the extracted value
may default to that placeholder.

If a pattern contains multiple captures, the recipe should specify which one is
the result, for example:

```yaml
extract:
  source: stdout
  match:
    mode: line
    pattern: "Energy = %(float:energy) Variance = %(float:variance)"
  value: energy
```

## Built-In Extractors

To avoid exposing even the pattern language in simple cases, the runtime should
ship with small typed extractors such as:

- `last_float`
- `first_float`
- `single_float`
- `json_field`
- `key_value`

That allows recipes like:

```yaml
recipes:
  norm:
    probe: mp-norm {state}
    extract: last_float
```

## Regex as Escape Hatch

Raw regex should still be available for awkward tools, but it should not be the
main authoring interface.

If regex is needed, it should live in the recipe library, not in every test.
