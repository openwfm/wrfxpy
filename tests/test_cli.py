from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def load_cli():
    loader = SourceFileLoader("wrfx_cli", str(ROOT / "wrfx"))
    spec = spec_from_loader(loader.name, loader)
    module = module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_primary_commands_dispatch_to_existing_entry_points(monkeypatch):
    cli = load_cli()
    calls = []
    monkeypatch.setattr(cli, "run_script", lambda script, arguments: calls.append((script, arguments)) or 0)

    assert cli.main(["forecast", "jobs/example.json"]) == 0
    assert calls.pop() == ("forecast.py", [str((cli.CALLER_DIR / "jobs/example.json").resolve())])

    assert cli.main(["list"]) == 0
    assert calls.pop() == ("cleanup.py", ["list"])

    assert cli.main(["cancel", "job-1"]) == 0
    assert calls.pop() == ("cleanup.py", ["cancel", "job-1"])

    assert cli.main(["clean", "output", "job-1"]) == 0
    assert calls.pop() == ("cleanup.py", ["output", "job-1"])


def test_cfsr_uses_the_generic_grib_command(monkeypatch):
    cli = load_cli()
    calls = []
    monkeypatch.setattr(cli, "run_script", lambda script, arguments: calls.append((script, arguments)) or 0)

    assert cli.main([
        "ingest", "gribs", "CFSR_P",
        "2016-09-10_00:00:00", "2016-09-13_00:00:00", "ingest",
    ]) == 0
    assert calls == [(
        "ingest/retrieve_gribs.py",
        [
            "CFSR_P",
            "2016-09-10_00:00:00",
            "2016-09-13_00:00:00",
            str((cli.CALLER_DIR / "ingest").resolve()),
        ],
    )]


def test_specialized_command_argument_order(monkeypatch):
    cli = load_cli()
    calls = []
    monkeypatch.setattr(cli, "run_script", lambda script, arguments: calls.append((script, arguments)) or 0)

    assert cli.main(["grib", "to-netcdf", "input.grib2", "7", "output.nc"]) == 0
    assert calls.pop() == (
        "ingest/grib_file.py",
        [
            "to_netcdf",
            str((cli.CALLER_DIR / "input.grib2").resolve()),
            "7",
            str((cli.CALLER_DIR / "output.nc").resolve()),
        ],
    )

    assert cli.main(["ingest", "landfire", "-112.8", "-112.1", "39.4", "39.9"]) == 0
    assert calls.pop() == (
        "ingest/retrieve_landfire.py",
        ["-112.8,-112.1,39.4,39.9"],
    )

    assert cli.main(["process", "file", "wrfout", "T2,PSFC", "products/frame", "2"]) == 0
    assert calls.pop() == (
        "vis/postprocess.py",
        [
            str((cli.CALLER_DIR / "wrfout").resolve()),
            "T2,PSFC",
            str((cli.CALLER_DIR / "products/frame").resolve()),
            "2",
        ],
    )

    assert cli.main([
        "geogrid", "convert", "fuel.tif", "geo_data", "NFUEL_CAT",
        "--bounds", "-112.8", "-112.1", "39.4", "39.9",
    ]) == 0
    assert calls.pop() == (
        "geo/convert_geotiff.py",
        [
            str((cli.CALLER_DIR / "fuel.tif").resolve()),
            str((cli.CALLER_DIR / "geo_data").resolve()),
            "NFUEL_CAT",
            "-112.8,-112.1,39.4,39.9",
        ],
    )

    assert cli.main(["process", "file", "wrfout", "@variables.json", "products/frame"]) == 0
    assert calls.pop()[1][1] == "@" + str((cli.CALLER_DIR / "variables.json").resolve())


def test_cycle_coordinate_forms_parse():
    cli = load_cli()
    parser = cli.build_parser()
    hrrr = parser.parse_args(["cycle", "hrrr", "a", "42", "-124.6", "49", "-116.4"])
    assert hrrr.arguments == ["a", "42", "-124.6", "49", "-116.4"]
    rtma = parser.parse_args(["cycle", "rtma", "42", "-124.6", "49", "-116.4"])
    assert rtma.arguments == ["42", "-124.6", "49", "-116.4"]


def test_arcgis_uses_its_existing_conda_environment(monkeypatch):
    cli = load_cli()
    calls = []
    monkeypatch.setattr(
        cli.subprocess,
        "call",
        lambda command, cwd, env: calls.append((command, cwd, env)) or 0,
    )

    assert cli.main(["ingest", "arcgis", "-112.8", "-112.1", "39.4", "39.9", "/tmp/fire"]) == 0
    command, cwd, environment = calls[0]
    assert command[:5] == ["conda", "run", "--no-capture-output", "-n", "arcgis"]
    assert command[-5:] == ["-112.8", "-112.1", "39.4", "39.9", "/tmp/fire"]
    assert cwd == str(cli.INSTALL_DIR)
    assert environment["PYTHONPATH"] == str(cli.SOURCE_DIR)


def test_satellite_processing_preserves_the_two_step_workflow(monkeypatch):
    cli = load_cli()
    calls = []
    monkeypatch.setattr(cli, "run_script", lambda script, arguments: calls.append((script, arguments)) or 0)

    assert cli.main(["process", "satellite", "job-1"]) == 0
    assert calls == [
        ("process_sat_output.py", ["job-1"]),
        ("make_kmz.py", ["job-1"]),
    ]


def test_csv_to_kmz_preserves_doc_kml_archive_name(monkeypatch, tmp_path):
    cli = load_cli()

    def create_kml(script, arguments):
        assert script == "vis/csv2kml.py"
        Path(arguments[1]).write_text("<kml/>")
        return 0

    monkeypatch.setattr(cli, "run_script", create_kml)
    csv_path = tmp_path / "points.csv"
    output_path = tmp_path / "points.kmz"
    csv_path.write_text("name,lat,lon\n")

    assert cli.main(["kmz", "from-csv", str(csv_path), str(output_path)]) == 0
    with zipfile.ZipFile(str(output_path)) as archive:
        assert archive.namelist() == ["doc.kml"]
        assert archive.read("doc.kml") == b"<kml/>"


def test_help_does_not_import_scientific_modules(capsys):
    cli = load_cli()
    try:
        cli.main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0
    assert "wrfx forecast" not in capsys.readouterr().err
