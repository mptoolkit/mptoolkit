# Contributed Models

Models in this directory are contributed, experimental, specialized, or otherwise
not part of the supported model set. They are useful examples and may be used in
research workflows, but they are not maintained to the same stability standard
as the generators in `models/`.

Conventions:

- `make contrib-models` builds every `*.cpp` file directly in this directory.
- `make all-contrib` builds the supported models, contributed models, and tools.
- Files under `models/contrib/non-functional/` are intentionally excluded from
  auto-discovery and may be incomplete or kept only for reference.
- Contributed models are not required to have integration tests.
- Promotion from `models/contrib/` to `models/` should include at least
  contract/smoke coverage and a convention check for descriptions, operators,
  and help text.
- Physics-sign or Hamiltonian-convention changes should not be made casually;
  add a small reference calculation first where practical.
