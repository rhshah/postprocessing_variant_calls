import pytest  # type: ignore
import os 
from typer.testing import CliRunner
from pdb import set_trace as bp
from postprocessing_variant_calls.main import app

runner = CliRunner()
vardict_single_calls = [
        ['vardict', 'single', 'filter', '--inputVcf', 'tests/data/vardict/single_test.vcf' , '--tsampleName', 'Myeloid200-1', '-ad','1', '-o', 'tests/data/vardict/single'], 

]

vardict_matched = [
    ['vardict', 'case-control', 'filter', '--inputVcf', 'tests/data/vardict/case_control_test.vcf' , '--tsampleName', 'C-C1V52M-L001-d', '-ad','1' , '-o', 'tests/data/vardict/two']
]

@pytest.mark.parametrize("call", vardict_single_calls)
def test_single(call):
    # vardict_tests['single']
    result = runner.invoke(app, call)
    result.stdout 
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/data/vardict/single/single_test_complex_STDfilter.vcf") == True
    assert os.path.exists("tests/data/vardict/single/single_test_STDfilter.vcf") == True
    assert os.path.exists("tests/data/vardict/single/single_test_STDfilter.txt") == True
    os.remove("tests/data/vardict/single/single_test_complex_STDfilter.vcf")
    os.remove("tests/data/vardict/single/single_test_STDfilter.vcf")
    os.remove("tests/data/vardict/single/single_test_STDfilter.txt")

@pytest.mark.parametrize("call", vardict_matched)
def test_two(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/data/vardict/two/case_control_test_complex_STDfilter.vcf") == True
    assert os.path.exists("tests/data/vardict/two/case_control_test_STDfilter.vcf") == True
    assert os.path.exists("tests/data/vardict/two/case_control_test_STDfilter.txt") == True
    os.remove("tests/data/vardict/two/case_control_test_complex_STDfilter.vcf")
    os.remove("tests/data/vardict/two/case_control_test_STDfilter.vcf")
    os.remove("tests/data/vardict/two/case_control_test_STDfilter.txt")


