# Suggested command-line improvements

The `wrfx` interface intentionally preserves the former shell script names,
arguments, working-directory behavior, and command sequences. The following
changes may improve the interface, but they are not part of the compatibility
refactoring:

- Group related commands, such as `wrfx ingest gribs` and
  `wrfx process output`, after users have an explicit migration plan.
- Resolve user-supplied paths relative to the invocation directory
  consistently. The former scripts do not all use the same rule.
- Replace comma-separated geographic bounds with four validated numeric
  arguments.
- Expose the optional postprocessing skip argument. `postprocess.sh` accepted
  only its first three arguments even though the Python entry point supports a
  fourth argument.
- Replace the external `zip` command in `csv2kmz` with Python's `zipfile`
  module and avoid leaving `doc.kml` in the installation directory.
- Replace shell-based ArcGIS environment activation with a documented,
  directly executable environment command.
- Convert the Python command entry points to importable `main(argv)` functions
  instead of executing their `__main__` blocks through `runpy`.
- Add consistent argument validation and per-command help.
