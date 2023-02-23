import pytest  # type: ignore
import os 
from typer.testing import CliRunner
from pdb import set_trace as bp
from postprocessing_variant_calls.main import app

runner = CliRunner()
maf_concat_files = [
        ['maf', 'concat', '-f', 'tests/data/maf/concat/C-1234-L001-d.maf' , 
        '-f', 'tests/data/maf/concat/C-1234-L002-d.maf', 
        "-o", "tests/data/maf/concat/output_maf"] 

]

maf_concat_paths = [
        ['maf', 'concat', '-p', 'tests/data/maf/concat/paths.txt' ,  
        "-o", "tests/data/maf/concat/output_maf"] 

]

@pytest.mark.parametrize("call", maf_concat_files)
def test_concat_files(call):
    result = runner.invoke(app, call)
    result.stdout 
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/data/maf/concat/output_maf.maf") == True
    os.remove("tests/data/maf/concat/output_maf.maf")

@pytest.mark.parametrize("call", maf_concat_paths)
def test_concat_paths(call):
    result = runner.invoke(app, call)
    result.stdout 
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/data/maf/concat/output_maf.maf") == True
    os.remove("tests/data/maf/concat/output_maf.maf")