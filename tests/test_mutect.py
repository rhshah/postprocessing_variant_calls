import pytest  # type: ignore
import os
from typer.testing import CliRunner
from pdb import set_trace as bp
from postprocessing_variant_calls.main import app

runner = CliRunner()

mutect1_filter_calls = [["mutect1", "case-control", "filter", "--help"]]


# Command Test
@pytest.mark.parametrize("call", mutect1_filter_calls)
def test_mutect1_filter(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0
