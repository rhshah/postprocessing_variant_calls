#!/usr/bin/env python
import pytest  # type: ignore
import os 
from typer.testing import CliRunner

from pv import app

runner = CliRunner()
vardict_single_calls = [
        ['vardict', 'single', 'filter', '--inputVcf', 'data/Myeloid200-1.vcf' , '--tsampleName', 'Myeloid200-1', '-ad','1', '-o', 'data/single'], 

]

vardict_matched = [
    ['vardict', 'case-control', 'filter', '--inputVcf', 'data/C-C1V52M-L001-d.DONOR22-TP.vardict.vcf' , '--tsampleName', 'C-C1V52M-L001-d', '-ad','1' , '-o', 'data/two']
]

@pytest.mark.parametrize("call", vardict_single_calls)
def test_single(call):
    # vardict_tests['single']
    result = runner.invoke(app, call)
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("data/single/Myeloid200-1_complex_STDfilter.vcf") == True
    assert os.path.exists("data/single/Myeloid200-1_STDfilter.vcf") == True
    assert os.path.exists("data/single/Myeloid200-1_STDfilter.txt") == True
    os.remove("data/single/Myeloid200-1_complex_STDfilter.vcf")
    os.remove("data/single/Myeloid200-1_STDfilter.vcf")
    os.remove("data/single/Myeloid200-1_STDfilter.txt")

@pytest.mark.parametrize("call", vardict_single_calls)
def test_two(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_complex_STDfilter.vcf") == True
    assert os.path.exists("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_STDfilter.vcf") == True
    assert os.path.exists("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_STDfilter.txt") == True
    os.remove("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_complex_STDfilter.vcf")
    os.remove("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_STDfilter.vcf")
    os.remove("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_STDfilter.txt")


