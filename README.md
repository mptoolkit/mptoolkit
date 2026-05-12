# Matrix Product Toolkit

This is the Github repository for the Matrix Product Toolkit.

The Matrix Product Toolkit is a project which started life around 2002, originally envisaged as a ‘next generation’ DMRG code, incorporating non-abelian symmetries and with an emphasis on a flexible and generic way to construct Hamiltonian operators and measure observables.  It has been used in over 100 research publications.

Documentation for the toolkit is available from https://mptoolkit.qusim.net/

## CMake build bootstrap

MPToolkit includes an initial CMake build. It is intended to live alongside the
existing Autoconf build while the CMake configuration is being completed.

Configure an out-of-tree build with a user-writable install prefix:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/mptk"
```

The CMake build defaults to `Release`. For a developer debug build, configure a
separate build tree:

```sh
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
```

Build and install the standard tools and model generators:

```sh
cmake --build build --target mptoolkit --parallel
cmake --install build --component tools
cmake --install build --component models
```

With the prefix above, executables are installed under `$HOME/mptk/bin`. Add
that directory to your `PATH`, for example:

```sh
export PATH="$HOME/mptk/bin:$PATH"
```

Other common non-root prefixes are `$HOME/.local`, which installs executables
under `$HOME/.local/bin`, and `$HOME`, which installs executables under
`$HOME/bin`.
