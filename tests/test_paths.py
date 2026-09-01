"""
Tests for the dataset resolver.

None of the source datasets ship with this code, so the failure a new user meets first is a
missing input. It has to name the dataset and the path it expected, not fail somewhere inside a
read. paths.py resolves ARAVALLI_WA_DATA at import time, so each test reloads it under a
controlled environment.

Author: Bhavik Harish Lodhia, Curtin University
"""
import importlib
import os

import pytest


def paths_rooted_at(directory, monkeypatch):
    monkeypatch.setenv("ARAVALLI_WA_DATA", str(directory))
    import paths
    return importlib.reload(paths)


class TestRequire:
    def test_raises_when_the_data_directory_is_empty(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        with pytest.raises(SystemExit):
            paths.require("NGSA")

    def test_message_names_the_dataset_and_the_expected_path(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        with pytest.raises(SystemExit) as excinfo:
            paths.require("NGSA")
        message = str(excinfo.value)
        assert "NGSA" in message
        assert "NGSA_data.csv" in message
        assert "data/README.md" in message

    def test_message_lists_every_missing_dataset_not_just_the_first(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        with pytest.raises(SystemExit) as excinfo:
            paths.require("NGSA", "HMMA", "NGCM_TABLE")
        message = str(excinfo.value)
        assert "NGSA" in message and "HMMA" in message and "NGCM_TABLE" in message

    def test_passes_silently_when_the_file_is_present(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        os.makedirs(os.path.dirname(paths.NGSA), exist_ok=True)
        open(paths.NGSA, "w").close()
        paths.require("NGSA")            # must not raise

    def test_an_unknown_name_is_ignored_rather_than_crashing(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        paths.require("NOT_A_DATASET")   # must not raise

    def test_no_arguments_is_a_no_op(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        paths.require()


class TestResolution:
    def test_the_environment_variable_overrides_the_repository_data_directory(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        assert paths.DATA == str(tmp_path)
        assert paths.NGSA.startswith(str(tmp_path))

    def test_results_stay_inside_the_repository_when_data_is_moved_out(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        assert not paths.RESULTS.startswith(str(tmp_path))
        assert paths.RESULTS.endswith("results")
        assert paths.DRAINAGE.endswith(os.path.join("results", "drainage"))

    def test_every_required_dataset_has_a_resolved_path(self, tmp_path, monkeypatch):
        paths = paths_rooted_at(tmp_path, monkeypatch)
        for name, path in paths._REQUIRED.items():
            assert isinstance(path, str) and path, name
