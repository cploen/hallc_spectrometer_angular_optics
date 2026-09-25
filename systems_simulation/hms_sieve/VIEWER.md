# Layered viewer: geometry, names and diagnostics

[Overview](README.md) · [iFarm build guide](IFARM.md) · [Reported baseline validation](IFARM_VALIDATION.md)

One immutable scene supplies every view. This update changes annotations and
visibility only: geometry parameters, transforms, equations, the A/B/C
constructions and all generated scene/numerical files remain unchanged.

## Start a clean view on iFarm

From the repository root, after obtaining this revision:

```sh
module load geant4/11.2.1
hms_repo_dir=$PWD
hms_build_dir=$(mktemp -d)
cmake -S "$hms_repo_dir/systems_simulation/hms_sieve" -B "$hms_build_dir" \
  -DHMS_WITH_INTERACTIVE=ON
cmake --build "$hms_build_dir" --parallel 2
cd "$hms_build_dir"
./hms_geometry_viewer \
  "$hms_repo_dir/systems_simulation/hms_sieve/generated/scene.tsv" \
  view.mac --interactive
```

Run from the **build directory**: CMake copies all macros there, and presets
call `reset_layers.mac` by its relative name. Reconfigure CMake when updating an
older build so the new macros are copied. No geometry regeneration is needed.

The default shows the beam/HMS axes, target and foil markers, reference planes,
and nominal sieve grid. It has only short plane/axis names. Triads, naming keys
and A/B/C diagnostics start hidden.

## Choose a view

In the Geant4 **Session** command box, enter, for example:

```text
/control/execute view_lab_frame.mac
```

These presets operate on the already-open viewer. They reset visibility, then
enable the named layers; they preserve your camera/orbit. Start with `view.mac`,
not a preset, when launching a new process.

| Preset macro | Visible content in addition to base geometry |
|---|---|
| `view_geometry.mac` | Short plane/axis names; clean default |
| `view_lab_frame.mac` | LAB triad |
| `view_hms_frame.mac` | HMS TRANSPORT triad |
| `view_sieve_frame.mac` | Sieve-local triad |
| `view_frames.mac` | All three triads |
| `view_hcana.mac` | Code naming key, including constructed quantities |
| `view_branches.mac` | Direct code-to-ROOT aliases |
| `view_aliases.mac` | Code and ROOT naming keys side by side |
| `view_rays.mac` | A/B/C rays, intersection markers, residual connectors |
| `view_everything.mac` | All layers; useful for inventory, not the clean overview |

Every preset includes the short plane/axis names. Hide those separately if
desired. Individual named Geant4 models can be switched without resetting the
other layers:

```text
/vis/scene/activateModel HMS_lab_frame true
/vis/scene/activateModel HMS_transport_frame false
/vis/scene/activateModel HMS_sieve_frame false
/vis/scene/activateModel HMS_labels false
/vis/viewer/refresh
```

| Model name | Content |
|---|---|
| `HMS_base` | Existing geometry primitives; short study/status footer |
| `HMS_labels` | Plane names and beam/HMS axis names |
| `HMS_lab_frame` | LAB axes and short labels |
| `HMS_transport_frame` | HMS axes and short labels |
| `HMS_sieve_frame` | Sieve-local axes and short labels |
| `HMS_code_names` | Code names and meaning; separately identified constructed names |
| `HMS_branch_names` | ROOT alias key |
| `HMS_rays` | Existing diagnostic rays and sieve intersections, short A/B/C labels |

Use `/vis/scene/list` to inspect model names. The Qt scene tree can also expose
these named models. Avoid `/vis/scene/activateModel all false`: Geant4 rejects
that command. Use the specific model names or `view_geometry.mac` instead.

## How to read the names

The naming keys stay in screen space instead of accumulating long labels at
nearby 3D endpoints. Each alias row is generated from one shared quantity record:
for example, `xtar = H.gtr.x`. **`=` means the same stored quantity**, not agreement
between different ray constructions. `[T]` identifies the HMS target plane and
direction convention; `[S]` identifies the sieve projection plane. Those tags
appear on the same geometric concepts in both naming views.

`[V]` means lab reaction coordinates; `[R]` means raster/BPM quantities. Their
provider mapping/reference remains OPEN. No saved reaction point or BPM plane
is fabricated just to attach a label. The scene contains **no saved event
values**: ROOT branches are a naming key, not assignments to the synthetic rays.

Keep these distinctions:

- `xtar`, `ytar`, `xptar`, `yptar` are saved names; their `T`-suffixed counterparts
  are the quantities constructed for fitting in B.
- `ztar = reactz = H.react.z`; `ztarT` is a nominal foil coordinate.
- `xsieve/ysieve` name saved projections; `xsT/ysT` are assigned nominal centers.
- `->` denotes construction from inputs, not an alias.
- Lab `Z=0` and HMS `z=0` are different planes. The sieve projection plane is
  defined; physical plate faces remain unresolved.

The diagnostic view labels A as the endpoint ray, B as the fit construction,
and C as the HCANA coordinate step with analytic inputs. Gray connectors join
the existing A/B/C sieve markers; they do not shift the markers or modify rays.
C is not a replay or a stored ROOT event. Numeric residuals remain in the
unchanged `generated/closure.tsv`.

## Verification status and short native check

The earlier single-layer viewer was built and displayed on iFarm at `7fff412`,
as reported in the handoff. That validates the baseline build path. **This
layered revision still needs its own iFarm build/display check.** Locally, the
11 existing geometry tests and four annotation checks pass. The new checks
compile the shared C++ annotation metadata, verify all 114 primitives are
assigned exactly once, check aliases and preset switching, and confirm that
geometry source/configuration hashes match the frozen report. They do not
substitute for Geant4 rendering.

After rebuilding on iFarm:

1. Confirm the default view has no triads, long variable names or A/B/C rays.
2. Switch each single-frame preset; only the selected triad should appear.
3. Compare the code and branch keys; confirm `xtarT` is not presented as a ROOT
   alias and that no key labels C as a saved event.
4. Show/hide `HMS_rays`; inspect the existing intersections without changing
   the underlying grid. Return to `view_geometry.mac`; the camera should stay put.
5. Save clean overview, single-frame and alias-view screenshots with the tested
   commit and Geant4 version. Record any label overlap or unsupported command.

`export_vrml.mac` exports the geometry, frames and diagnostics. The iFarm
`VRML2FILE` driver does not support `G4Text`; inspect names in Qt/OpenGL, not in
the VRML output. That known driver limitation is not a failed geometry export.

Implementation uses Geant4's named user actions and
[`/vis/scene/activateModel`](https://raw.githubusercontent.com/Geant4/geant4/master/source/visualization/management/src/G4VisCommandsScene.cc),
with [user-action descriptions as model names](https://raw.githubusercontent.com/Geant4/geant4/master/source/visualization/management/src/G4VisCommandsSceneAdd.cc).
No custom GUI framework or second geometry is introduced.
