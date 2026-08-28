from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_cli():
    loader = SourceFileLoader("wrfx_cli", str(ROOT / "wrfx"))
    spec = spec_from_loader(loader.name, loader)
    module = module_from_spec(spec)
    loader.exec_module(module)
    return module


def capture_scripts(monkeypatch, cli):
    calls = []
    monkeypatch.setattr(
        cli,
        "run_script",
        lambda script, arguments, work_dir=None: calls.append((script, arguments, work_dir)) or 0,
    )
    return calls


def test_commands_keep_shell_script_basenames():
    cli = load_cli()
    assert cli.COMMANDS == (
        "apply_fmda", "clamp2mesh", "cleanup", "convert_geotiff", "csv2kmz",
        "domain_bbox", "domain_setup", "earthdata", "execute_wrf", "fill_subgrid",
        "forecast", "grib_tool", "hrrr_cycler", "level0_retr", "make_job_file",
        "make_kmz", "postprocess", "process_output", "process_sat_output",
        "process_tiffs_output", "recover_catalog", "retrieve_arcgis", "retrieve_gribs",
        "retrieve_hdfs", "retrieve_landfire", "retrieve_sat", "rtma_cycler",
        "simple_forecast", "ssh_command", "ssh_shuttle",
    )


def test_primary_commands_keep_arguments_and_working_directory(monkeypatch):
    cli = load_cli()
    calls = capture_scripts(monkeypatch, cli)

    assert cli.main(["forecast", "jobs/example.json", "ignored"]) == 0
    assert calls.pop() == ("forecast.py", ["jobs/example.json"], cli.INSTALL_DIR)

    assert cli.main(["cleanup", "list"]) == 0
    assert calls.pop() == ("cleanup.py", ["list"], cli.INSTALL_DIR)

    assert cli.main(["cleanup", "output", "job-1", "ignored"]) == 0
    assert calls.pop() == ("cleanup.py", ["output", "job-1"], cli.INSTALL_DIR)

    assert cli.main(["process_output", "job-1", "ignored"]) == 0
    assert calls.pop() == ("process_output.py", ["job-1"], cli.INSTALL_DIR)


def test_retrieval_commands_keep_argument_forms(monkeypatch):
    cli = load_cli()
    calls = capture_scripts(monkeypatch, cli)

    grib_arguments = ["CFSR_P", "2016-09-10_00:00:00", "2016-09-13_00:00:00", "ingest"]
    assert cli.main(["retrieve_gribs"] + grib_arguments + ["ignored"]) == 0
    assert calls.pop() == ("ingest/retrieve_gribs.py", grib_arguments, cli.INSTALL_DIR)

    landfire_arguments = ["-112.8,-112.1,39.4,39.9"]
    assert cli.main(["retrieve_landfire"] + landfire_arguments) == 0
    assert calls.pop() == ("ingest/retrieve_landfire.py", landfire_arguments, cli.INSTALL_DIR)

    hdfs_arguments = [
        "MODIS_AQUA", "2017.06.20-10.00.00", "2017.06.20-15.00.00", "hdf_test",
        "-124.7", "-66.9", "24.7", "49.3",
    ]
    assert cli.main(["retrieve_hdfs"] + hdfs_arguments) == 0
    assert calls.pop() == ("ingest/retrieve_hdfs.py", hdfs_arguments, None)


def test_special_shell_argument_behavior_is_preserved(monkeypatch):
    cli = load_cli()
    calls = capture_scripts(monkeypatch, cli)

    assert cli.main(["make_kmz", "job", "extra"]) == 0
    assert calls.pop() == ("make_kmz.py", ["job extra"], cli.INSTALL_DIR)

    assert cli.main(["postprocess", "wrfout", "T2,PSFC", "prefix", "2"]) == 0
    assert calls.pop() == ("vis/postprocess.py", ["wrfout", "T2,PSFC", "prefix"], None)

    assert cli.main(["process_sat_output", "job-1"]) == 0
    assert calls == [
        ("process_sat_output.py", ["job-1"], cli.INSTALL_DIR),
        ("make_kmz.py", ["job-1"], cli.INSTALL_DIR),
    ]


def test_convert_geotiff_keeps_caller_relative_path_conversion(monkeypatch, tmp_path):
    cli = load_cli()
    calls = capture_scripts(monkeypatch, cli)
    monkeypatch.chdir(tmp_path)

    arguments = ["fuel.tif", "geo_data", "NFUEL_CAT", "-112.8,-112.1,39.4,39.9"]
    assert cli.main(["convert_geotiff"] + arguments) == 0
    assert calls == [(
        "geo/convert_geotiff.py",
        [str(tmp_path / "fuel.tif"), str(tmp_path / "geo_data"), "NFUEL_CAT", arguments[3]],
        cli.INSTALL_DIR,
    )]


def test_csv2kmz_keeps_external_zip_workflow(monkeypatch):
    cli = load_cli()
    calls = capture_scripts(monkeypatch, cli)
    zip_calls = []
    monkeypatch.setattr(cli.os, "remove", lambda path: None)
    monkeypatch.setattr(
        cli.subprocess,
        "call",
        lambda command, cwd=None: zip_calls.append((command, cwd)) or 0,
    )

    assert cli.main(["csv2kmz", "points.csv", "points.kmz"]) == 0
    assert calls == [("vis/csv2kml.py", ["points.csv", "doc.kml"], cli.INSTALL_DIR)]
    assert zip_calls == [(["zip", "-9", "points.kmz", "doc.kml"], str(cli.INSTALL_DIR))]


def test_arcgis_keeps_bashrc_and_conda_activation(monkeypatch):
    cli = load_cli()
    calls = []
    monkeypatch.setattr(
        cli.subprocess,
        "call",
        lambda command, cwd=None: calls.append((command, cwd)) or 0,
    )
    arguments = ["-112.8", "-112.1", "39.4", "39.9", "/tmp/fire"]

    assert cli.main(["retrieve_arcgis"] + arguments) == 0
    command, cwd = calls[0]
    assert command[:3] == ["bash", "-lc", command[2]]
    assert "source ~/.bashrc; conda activate arcgis" in command[2]
    assert command[-5:] == arguments
    assert cwd is None


def test_commands_that_print_pwd_keep_doing_so(monkeypatch, capsys):
    cli = load_cli()
    calls = capture_scripts(monkeypatch, cli)

    assert cli.main(["hrrr_cycler", "a", "FIRE"]) == 0
    assert calls.pop() == ("hrrr_cycler.py", ["a", "FIRE"], None)
    assert str(Path.cwd()) in capsys.readouterr().out


def test_missing_required_argument_uses_command_name():
    cli = load_cli()
    try:
        cli.main(["forecast"])
    except SystemExit as exc:
        assert str(exc) == "usage: wrfx forecast input.json"
    else:
        raise AssertionError("forecast without an input should fail")
