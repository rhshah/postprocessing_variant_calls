import pytest  # type: ignore
import os 
from typer.testing import CliRunner
from pdb import set_trace as bp
from postprocessing_variant_calls.main import app

runner = CliRunner()
maf_concat_files = [
        ['maf', 'concat', '-f', 'tests/data/maf/concat/maf2.maf' , 
        '-f', 'tests/data/maf/concat/maf1.maf', 
        "-o", "tests/data/maf/concat/output_maf.maf", "-h" ,"resources/maf_concat/header.txt"] 

]

maf_concat_paths = [
        ['maf', 'concat', '-p', 'tests/data/maf/concat/paths.txt' ,  
        "-o", "tests/data/maf/concat/output_maf.maf", "-h", "resources/maf_concat/header.txt"] 

]

maf_annotate_maf_by_bed = [
        ['maf', 'annotate', 'mafbybed', "--help"] 

]

maf_annotate_maf_by_bed = [
        ['maf', 'annotate', 'mafbybed', "--help"] 
        ]

maf_subset = [
    ['maf', 'subset', '-m', 'tests/data/maf/subset/example_input.maf', '--ids', 'resources/maf_subset/example_subset_ids.txt', '-o','tests/data/maf/subset/output_subset.maf']
    ]

maf_annotate_help_message = 'Usage: root maf annotate mafbybed [OPTIONS]\n\nOptions:\n  -m, --maf FILE     input maf file  [required]\n  -b, --bed FILE     bed file to annotate maf  [required]\n  -o, --output TEXT  output maf file  [default: output]\n  -c, --cname TEXT   name for annotation column  [default: annotation]\n  --help             Show this message and exit.\n'

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


@pytest.mark.parametrize("call", maf_annotate_maf_by_bed)
def test_annotate_mafbybed(call):
    result = runner.invoke(app, call)
    assert result.stdout == maf_annotate_help_message


@pytest.mark.parametrize("call", maf_subset)
def test_concat_files(call):
    result = runner.invoke(app, call)
    result.stdout 
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/data/maf/subset/output_subset.maf") == True
    os.remove("tests/data/maf/subset/output_subset.maf")
