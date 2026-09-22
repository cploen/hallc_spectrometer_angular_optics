# Using the geometry study on iFarm

[Overview and file guide](README.md) · [Testing](TESTING.md) · [Next steps](NEXT_STEPS.md)

These commands use the full repository layout. **They have not been executed
on iFarm, and the native viewer has not yet compiled/run locally.** No particular
iFarm module name, Geant4 installation path or display service is assumed.

## Transfer the committed work

A local commit is not automatically available on GitHub. Once the branch has
been pushed, a new checkout can use:

```sh
git clone --branch optics-validation-studies \
  https://github.com/cploen/hallc_spectrometer_angular_optics.git
cd hallc_spectrometer_angular_optics
git log -1 --oneline
```

To transfer without pushing, create a bundle on the source machine:

```sh
git bundle create /tmp/hms-geometry.bundle optics-validation-studies
```

Copy that bundle using your usual transfer method. On iFarm, in a directory
where the new checkout should live, run:

```sh
git clone --branch optics-validation-studies /path/to/hms-geometry.bundle hallc-geometry
cd hallc-geometry
git log -1 --oneline
```

Replace the bundle path with its actual transferred location. Keep the full
checkout: the generator and parity test depend on parent repository files.
No external campaign ROOT files are needed.

## Validate coordinates first

Follow [TESTING.md](TESTING.md). These checks need only Python 3.7+ and a C++17
compiler; run them before configuring Geant4. Record the commit and the test
result, including any skipped tests.

## Configure and build Geant4 viewing support

Activate an existing Geant4 environment using the setup instructions for that
installation. Obtain the actual directory containing `Geant4Config.cmake`;
replace the placeholder in the command below. CMake 3.16+ and Geant4 visualization
libraries are required. A build without OpenGL can use the file renderer if its
Geant4 installation provides `VRML2FILE`.

From the repository root:

```sh
hms_repo_dir=$PWD
hms_build_dir=$(mktemp -d)
cmake -S "$hms_repo_dir/systems_simulation/hms_sieve" -B "$hms_build_dir" \
  -DGeant4_DIR=/actual/directory/containing/Geant4Config.cmake \
  -DHMS_WITH_INTERACTIVE=OFF
cmake --build "$hms_build_dir" --parallel 2
```

`Geant4_DIR` takes the **directory**, not the filename. If the installation's
setup already makes Geant4 discoverable to CMake, omit that option. Keep build
and generated output directories on a filesystem writable from the host where
you run these commands. No scheduler or site-specific environment setup is
invented here.

## File export and visual inspection

Use the committed scene first; run from the build directory so viewer outputs
do not clutter the repository:

```sh
cd "$hms_build_dir"
./hms_geometry_viewer \
  "$hms_repo_dir/systems_simulation/hms_sieve/generated/scene.tsv" \
  "$hms_build_dir/export_vrml.mac"
ls -l g4_*.wrl
```

Expected: exit status zero, a message identifying the loaded coordinate
primitives, and a nonempty `g4_*.wrl` file. Transfer the VRML file to a machine
with a compatible viewer if needed. Creating a file alone does not complete
visual validation: inspect the axes, two target planes, sieve grid and A/B/C
rays against the committed SVG and scene coordinates.

For a Geant4 build with interactive UI/OpenGL and a working graphical display:

```sh
cmake -S "$hms_repo_dir/systems_simulation/hms_sieve" -B "$hms_build_dir" \
  -DHMS_WITH_INTERACTIVE=ON
cmake --build "$hms_build_dir" --parallel 2
./hms_geometry_viewer \
  "$hms_repo_dir/systems_simulation/hms_sieve/generated/scene.tsv" \
  "$hms_build_dir/view.mac" --interactive
```

The same CMake cache retains `Geant4_DIR`. UI/OpenGL availability and remote
display access depend on the installed environment; file export is the first
check when no display is available. Neither mode runs events or studies
acceptance/scattering.

## Common stopping points

| Symptom | Check |
|---|---|
| CMake cannot find Geant4 | Activate the actual installation; locate its `Geant4Config.cmake` and pass the containing directory |
| Compiler parity test skips or cannot find the optics header | Use the full checkout and a C++17 compiler; set `CXX` to its executable if necessary |
| `/vis/open VRML2FILE` fails | Check that the installed Geant4 includes that visualization driver; record the version/build configuration |
| OpenGL/UI fails on a remote host | Use file mode first; resolve the site's display environment separately |
| Source hashes differ from the committed report | Confirm the commit, local edits and selected JSON configuration before interpreting residual changes |

Record the commit, configuration hash, host, compiler/CMake/Geant4 versions,
commands, output filenames, and visual-check result. Update the native-runtime
status only after actual success; see [NEXT_STEPS.md](NEXT_STEPS.md).
