#!/usr/bin/env python
import pytest  # type: ignore
import os 
from typer.testing import CliRunner

from postprocessing_variant_calls.src.vardict.vardict_process import app

runner = CliRunner()
sample = {
    'single':['--inputVcf','data/Myeloid200-1.vcf','--tsampleName','Myeloid200-1', '-ad', '1', '-o', 'data/single', '--single'], 
    'two':['--inputVcf','data/C-C1V52M-L001-d.DONOR22-TP.vardict.vcf','--tsampleName','C-C1V52M-L001-d', '-ad', '1', '-o', 'data/two', '--two']
}

def test_single():
    result = runner.invoke(app, sample['single'])
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("data/single/Myeloid200-1_complex_STDfilter.vcf") == True
    assert os.path.exists("data/single/Myeloid200-1_STDfilter.vcf") == True
    assert os.path.exists("data/single/Myeloid200-1_STDfilter.txt") == True

def test_two():
    result = runner.invoke(app, sample['two'])
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_complex_STDfilter.vcf") == True
    assert os.path.exists("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_STDfilter.vcf") == True
    assert os.path.exists("data/two/C-C1V52M-L001-d.DONOR22-TP.vardict_STDfilter.txt") == True


@pytest.mark.parametrize("name", sample)
def test_arg(name: str):
    # test single sample 
    test_single() 
    # test two sample
    test_two()

