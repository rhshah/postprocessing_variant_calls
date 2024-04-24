# Post-processing of variant calls

This package provides a variety of commands for manipulating different types of common outputs (e.g. mafs, vcf and txt files) from different bioinformatic variant callers such as mutect and vardict. 

Supported File Types:
- [maf](docs/MAF.md) 
- [vardict](docs/VARDICT.md)

# Installation  

For general use you can run: Collecting postprocessing_variant_calls
  Downloading postprocessing_variant_calls-0.2.8-py3-none-any.whl.metadata (2.0 kB)
Collecting PyVCF3 (from postprocessing_variant_calls)
  Downloading PyVCF3-1.0.3.tar.gz (977 kB)
     ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 977.6/977.6 kB 26.8 MB/s eta 0:00:00
  Installing build dependencies: started
  Installing build dependencies: finished with status 'done'
  Getting requirements to build wheel: started
  Getting requirements to build wheel: finished with status 'done'
  Installing backend dependencies: started
  Installing backend dependencies: finished with status 'done'
  Preparing metadata (pyproject.toml): started
  Preparing metadata (pyproject.toml): finished with status 'done'
Collecting numpy (from postprocessing_variant_calls)
  Downloading numpy-1.26.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (61 kB)
     ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 61.0/61.0 kB 18.1 MB/s eta 0:00:00
Collecting pandas (from postprocessing_variant_calls)
  Downloading pandas-2.2.2-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (19 kB)
Collecting typer[all] (from postprocessing_variant_calls)
  Downloading typer-0.12.3-py3-none-any.whl.metadata (15 kB)
Collecting python-dateutil>=2.8.2 (from pandas->postprocessing_variant_calls)
  Downloading python_dateutil-2.9.0.post0-py2.py3-none-any.whl.metadata (8.4 kB)
Collecting pytz>=2020.1 (from pandas->postprocessing_variant_calls)
  Downloading pytz-2024.1-py2.py3-none-any.whl.metadata (22 kB)
Collecting tzdata>=2022.7 (from pandas->postprocessing_variant_calls)
  Downloading tzdata-2024.1-py2.py3-none-any.whl.metadata (1.4 kB)
Requirement already satisfied: setuptools in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from PyVCF3->postprocessing_variant_calls) (65.5.0)
Collecting click>=8.0.0 (from typer[all]->postprocessing_variant_calls)
  Downloading click-8.1.7-py3-none-any.whl.metadata (3.0 kB)
Collecting typing-extensions>=3.7.4.3 (from typer[all]->postprocessing_variant_calls)
  Downloading typing_extensions-4.11.0-py3-none-any.whl.metadata (3.0 kB)
Requirement already satisfied: shellingham>=1.3.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from typer[all]->postprocessing_variant_calls) (1.5.4)
Collecting rich>=10.11.0 (from typer[all]->postprocessing_variant_calls)
  Downloading rich-13.7.1-py3-none-any.whl.metadata (18 kB)
Collecting six>=1.5 (from python-dateutil>=2.8.2->pandas->postprocessing_variant_calls)
  Downloading six-1.16.0-py2.py3-none-any.whl.metadata (1.8 kB)
Collecting markdown-it-py>=2.2.0 (from rich>=10.11.0->typer[all]->postprocessing_variant_calls)
  Downloading markdown_it_py-3.0.0-py3-none-any.whl.metadata (6.9 kB)
Collecting pygments<3.0.0,>=2.13.0 (from rich>=10.11.0->typer[all]->postprocessing_variant_calls)
  Downloading pygments-2.17.2-py3-none-any.whl.metadata (2.6 kB)
Collecting mdurl~=0.1 (from markdown-it-py>=2.2.0->rich>=10.11.0->typer[all]->postprocessing_variant_calls)
  Downloading mdurl-0.1.2-py3-none-any.whl.metadata (1.6 kB)
Downloading postprocessing_variant_calls-0.2.8-py3-none-any.whl (30 kB)
Downloading numpy-1.26.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (18.2 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 18.2/18.2 MB 58.1 MB/s eta 0:00:00
Downloading pandas-2.2.2-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (13.0 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 13.0/13.0 MB 69.6 MB/s eta 0:00:00
Downloading click-8.1.7-py3-none-any.whl (97 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 97.9/97.9 kB 29.7 MB/s eta 0:00:00
Downloading python_dateutil-2.9.0.post0-py2.py3-none-any.whl (229 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 229.9/229.9 kB 51.9 MB/s eta 0:00:00
Downloading pytz-2024.1-py2.py3-none-any.whl (505 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 505.5/505.5 kB 74.5 MB/s eta 0:00:00
Downloading rich-13.7.1-py3-none-any.whl (240 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 240.7/240.7 kB 51.4 MB/s eta 0:00:00
Downloading typing_extensions-4.11.0-py3-none-any.whl (34 kB)
Downloading tzdata-2024.1-py2.py3-none-any.whl (345 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 345.4/345.4 kB 45.5 MB/s eta 0:00:00
Downloading typer-0.12.3-py3-none-any.whl (47 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 47.2/47.2 kB 11.7 MB/s eta 0:00:00
Downloading markdown_it_py-3.0.0-py3-none-any.whl (87 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 87.5/87.5 kB 26.0 MB/s eta 0:00:00
Downloading pygments-2.17.2-py3-none-any.whl (1.2 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1.2/1.2 MB 66.9 MB/s eta 0:00:00
Downloading six-1.16.0-py2.py3-none-any.whl (11 kB)
Downloading mdurl-0.1.2-py3-none-any.whl (10.0 kB)
Building wheels for collected packages: PyVCF3
  Building wheel for PyVCF3 (pyproject.toml): started
  Building wheel for PyVCF3 (pyproject.toml): finished with status 'done'
  Created wheel for PyVCF3: filename=PyVCF3-1.0.3-py3-none-any.whl size=988536 sha256=f4070ddcde8f9edfd80e7ef6c1a565083cb74a049bd68d39f333771df8813760
  Stored in directory: /home/runner/.cache/pip/wheels/12/2e/72/c03483bc49d13b24ec16d0878e17df5e9d2cd9a088935acdda
Successfully built PyVCF3
Installing collected packages: pytz, tzdata, typing-extensions, six, PyVCF3, pygments, numpy, mdurl, click, python-dateutil, markdown-it-py, rich, pandas, typer, postprocessing_variant_calls
Successfully installed PyVCF3-1.0.3 click-8.1.7 markdown-it-py-3.0.0 mdurl-0.1.2 numpy-1.26.4 pandas-2.2.2 postprocessing_variant_calls-0.2.8 pygments-2.17.2 python-dateutil-2.9.0.post0 pytz-2024.1 rich-13.7.1 six-1.16.0 typer-0.12.3 typing-extensions-4.11.0 tzdata-2024.1
or a tagged version with 

For setting up a development environment please see the [Setting up a Dev Environment](#Setting-up-a-Dev-Environment) section.

# Usage

See [CLI](docs/CLI.md) for commmand line usage of the package.

# Setting up a Dev Environment 

## Install External Dependencies
Have an environment with python >= 3.8 installed. 

Install poetry: 

Requirement already satisfied: poetry in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (1.8.2)
Requirement already satisfied: build<2.0.0,>=1.0.3 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.2.1)
Requirement already satisfied: cachecontrol<0.15.0,>=0.14.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from cachecontrol[filecache]<0.15.0,>=0.14.0->poetry) (0.14.0)
Requirement already satisfied: cleo<3.0.0,>=2.1.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (2.1.0)
Requirement already satisfied: crashtest<0.5.0,>=0.4.1 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (0.4.1)
Requirement already satisfied: dulwich<0.22.0,>=0.21.2 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (0.21.7)
Requirement already satisfied: fastjsonschema<3.0.0,>=2.18.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (2.19.1)
Requirement already satisfied: installer<0.8.0,>=0.7.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (0.7.0)
Requirement already satisfied: keyring<25.0.0,>=24.0.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (24.3.1)
Requirement already satisfied: packaging>=23.1 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (24.0)
Requirement already satisfied: pexpect<5.0.0,>=4.7.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (4.9.0)
Requirement already satisfied: pkginfo<2.0.0,>=1.9.4 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.10.0)
Requirement already satisfied: platformdirs<5,>=3.0.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (4.2.1)
Requirement already satisfied: poetry-core==1.9.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.9.0)
Requirement already satisfied: poetry-plugin-export<2.0.0,>=1.6.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.7.1)
Requirement already satisfied: pyproject-hooks<2.0.0,>=1.0.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.0.0)
Requirement already satisfied: requests<3.0,>=2.26 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (2.31.0)
Requirement already satisfied: requests-toolbelt<2.0.0,>=1.0.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.0.0)
Requirement already satisfied: shellingham<2.0,>=1.5 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (1.5.4)
Requirement already satisfied: tomli<3.0.0,>=2.0.1 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (2.0.1)
Requirement already satisfied: tomlkit<1.0.0,>=0.11.4 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (0.12.4)
Requirement already satisfied: trove-classifiers>=2022.5.19 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (2024.4.10)
Requirement already satisfied: virtualenv<21.0.0,>=20.23.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from poetry) (20.26.0)
Requirement already satisfied: msgpack<2.0.0,>=0.5.2 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from cachecontrol<0.15.0,>=0.14.0->cachecontrol[filecache]<0.15.0,>=0.14.0->poetry) (1.0.8)
Requirement already satisfied: filelock>=3.8.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from cachecontrol[filecache]<0.15.0,>=0.14.0->poetry) (3.13.4)
Requirement already satisfied: rapidfuzz<4.0.0,>=3.0.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from cleo<3.0.0,>=2.1.0->poetry) (3.8.1)
Requirement already satisfied: urllib3>=1.25 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from dulwich<0.22.0,>=0.21.2->poetry) (2.2.1)
Requirement already satisfied: jaraco.classes in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from keyring<25.0.0,>=24.0.0->poetry) (3.4.0)
Requirement already satisfied: importlib-metadata>=4.11.4 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from keyring<25.0.0,>=24.0.0->poetry) (7.1.0)
Requirement already satisfied: SecretStorage>=3.2 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from keyring<25.0.0,>=24.0.0->poetry) (3.3.3)
Requirement already satisfied: jeepney>=0.4.2 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from keyring<25.0.0,>=24.0.0->poetry) (0.8.0)
Requirement already satisfied: ptyprocess>=0.5 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from pexpect<5.0.0,>=4.7.0->poetry) (0.7.0)
Requirement already satisfied: charset-normalizer<4,>=2 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from requests<3.0,>=2.26->poetry) (3.3.2)
Requirement already satisfied: idna<4,>=2.5 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from requests<3.0,>=2.26->poetry) (3.7)
Requirement already satisfied: certifi>=2017.4.17 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from requests<3.0,>=2.26->poetry) (2024.2.2)
Requirement already satisfied: distlib<1,>=0.3.7 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from virtualenv<21.0.0,>=20.23.0->poetry) (0.3.8)
Requirement already satisfied: zipp>=0.5 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from importlib-metadata>=4.11.4->keyring<25.0.0,>=24.0.0->poetry) (3.18.1)
Requirement already satisfied: cryptography>=2.0 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from SecretStorage>=3.2->keyring<25.0.0,>=24.0.0->poetry) (42.0.5)
Requirement already satisfied: more-itertools in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from jaraco.classes->keyring<25.0.0,>=24.0.0->poetry) (10.2.0)
Requirement already satisfied: cffi>=1.12 in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from cryptography>=2.0->SecretStorage>=3.2->keyring<25.0.0,>=24.0.0->poetry) (1.16.0)
Requirement already satisfied: pycparser in /opt/hostedtoolcache/Python/3.10.14/x64/lib/python3.10/site-packages (from cffi>=1.12->cryptography>=2.0->SecretStorage>=3.2->keyring<25.0.0,>=24.0.0->poetry) (2.22)

## Install Package Dependencies

Then install project dependencies with Poetry.



## Accessing Environment

To access the environment after initial setup up run: 



# Contributing to Documentation

The Gitbook for this repository is configured so changes are written in Gitbook and synced with the  branch. 

To contribute to the documentation, you can write your changes in [Gitbook](https://app.gitbook.com/o/-LhMNgvjydB3TFWAUMVb/s/VBp8SqbRAs28AQCVNoIS/), request a review, and merge the changes. Keep in mind, you will need access to the organization to contribute. 

Each file-type supported should have a section in the Gitbook detailing the implementation of the file-type and a justification of it's operations. For example, the  file-type has it's own section, which includes a description of how a maf is defined internally in the package and a justification of it's operations and how to use them.

Beyond file-type sections, you will also notice a section called Usage is: mono [options] program [program-options]

Development:
    --aot[=<options>]      Compiles the assembly to native code
    --debug[=<options>]    Enable debugging support, use --help-debug for details
    --debugger-agent=options Enable the debugger agent
    --profile[=profiler]   Runs in profiling mode with the specified profiler module
    --trace[=EXPR]         Enable tracing, use --help-trace for details
    --jitmap               Output a jit method map to /tmp/perf-PID.map
    --help-devel           Shows more options available to developers

Runtime:
    --config FILE          Loads FILE as the Mono config
    --verbose, -v          Increases the verbosity level
    --help, -h             Show usage information
    --version, -V          Show version information
    --version=number       Show version number
    --runtime=VERSION      Use the VERSION runtime, instead of autodetecting
    --optimize=OPT         Turns on or off a specific optimization
                           Use --list-opt to get a list of optimizations
    --attach=OPTIONS       Pass OPTIONS to the attach agent in the runtime.
                           Currently the only supported option is 'disable'.
    --llvm, --nollvm       Controls whenever the runtime uses LLVM to compile code.
    --gc=[sgen,boehm]      Select SGen or Boehm GC (runs mono or mono-sgen)
    --handlers             Install custom handlers, use --help-handlers for details.
    --aot-path=PATH        List of additional directories to search for AOT images., which lists all commands in the  package. Do not manually edit this section. This section is created using the  package, which uses the typer GNU bash, version 5.1.16(1)-release (x86_64-pc-linux-gnu)
These shell commands are defined internally.  Type `help' to see this list.
Type `help name' to find out more about the function `name'.
Use `info bash' to find out more about the shell in general.
Use `man -k' or `info' to find out more about commands not in this list.

A star (*) next to a name means that the command is disabled.

 job_spec [&]                            history [-c] [-d offset] [n] or hist>
 (( expression ))                        if COMMANDS; then COMMANDS; [ elif C>
 . filename [arguments]                  jobs [-lnprs] [jobspec ...] or jobs >
 :                                       kill [-s sigspec | -n signum | -sigs>
 [ arg... ]                              let arg [arg ...]
 [[ expression ]]                        local [option] name[=value] ...
 alias [-p] [name[=value] ... ]          logout [n]
 bg [job_spec ...]                       mapfile [-d delim] [-n count] [-O or>
 bind [-lpsvPSVX] [-m keymap] [-f file>  popd [-n] [+N | -N]
 break [n]                               printf [-v var] format [arguments]
 builtin [shell-builtin [arg ...]]       pushd [-n] [+N | -N | dir]
 caller [expr]                           pwd [-LP]
 case WORD in [PATTERN [| PATTERN]...)>  read [-ers] [-a array] [-d delim] [->
 cd [-L|[-P [-e]] [-@]] [dir]            readarray [-d delim] [-n count] [-O >
 command [-pVv] command [arg ...]        readonly [-aAf] [name[=value] ...] o>
 compgen [-abcdefgjksuv] [-o option] [>  return [n]
 complete [-abcdefgjksuv] [-pr] [-DEI]>  select NAME [in WORDS ... ;] do COMM>
 compopt [-o|+o option] [-DEI] [name .>  set [-abefhkmnptuvxBCHP] [-o option->
 continue [n]                            shift [n]
 coproc [NAME] command [redirections]    shopt [-pqsu] [-o] [optname ...]
 declare [-aAfFgiIlnrtux] [-p] [name[=>  source filename [arguments]
 dirs [-clpv] [+N] [-N]                  suspend [-f]
 disown [-h] [-ar] [jobspec ... | pid >  test [expr]
 echo [-neE] [arg ...]                   time [-p] pipeline
 enable [-a] [-dnps] [-f filename] [na>  times
 eval [arg ...]                          trap [-lp] [[arg] signal_spec ...]
 exec [-cl] [-a name] [command [argume>  true
 exit [n]                                type [-afptP] name [name ...]
 export [-fn] [name[=value] ...] or ex>  typeset [-aAfFgiIlnrtux] [-p] name[=>
 false                                   ulimit [-SHabcdefiklmnpqrstuvxPT] [l>
 fc [-e ename] [-lnr] [first] [last] o>  umask [-p] [-S] [mode]
 fg [job_spec]                           unalias [-a] name [name ...]
 for NAME [in WORDS ... ] ; do COMMAND>  unset [-f] [-v] [-n] [name ...]
 for (( exp1; exp2; exp3 )); do COMMAN>  until COMMANDS; do COMMANDS; done
 function name { COMMANDS ; } or name >  variables - Names and meanings of so>
 getopts optstring name [arg ...]        wait [-fn] [-p var] [id ...]
 hash [-lr] [-p pathname] [-dt] [name >  while COMMANDS; do COMMANDS; done
 help [-dms] [pattern ...]               { COMMANDS ; } parameters specified in typer commands to generate documentation. It is automatically updated by the git-action,  upon a push to the  branch. To make sure the  document updates to include newly added commands, specify all relevant typer GNU bash, version 5.1.16(1)-release (x86_64-pc-linux-gnu)
These shell commands are defined internally.  Type `help' to see this list.
Type `help name' to find out more about the function `name'.
Use `info bash' to find out more about the shell in general.
Use `man -k' or `info' to find out more about commands not in this list.

A star (*) next to a name means that the command is disabled.

 job_spec [&]                            history [-c] [-d offset] [n] or hist>
 (( expression ))                        if COMMANDS; then COMMANDS; [ elif C>
 . filename [arguments]                  jobs [-lnprs] [jobspec ...] or jobs >
 :                                       kill [-s sigspec | -n signum | -sigs>
 [ arg... ]                              let arg [arg ...]
 [[ expression ]]                        local [option] name[=value] ...
 alias [-p] [name[=value] ... ]          logout [n]
 bg [job_spec ...]                       mapfile [-d delim] [-n count] [-O or>
 bind [-lpsvPSVX] [-m keymap] [-f file>  popd [-n] [+N | -N]
 break [n]                               printf [-v var] format [arguments]
 builtin [shell-builtin [arg ...]]       pushd [-n] [+N | -N | dir]
 caller [expr]                           pwd [-LP]
 case WORD in [PATTERN [| PATTERN]...)>  read [-ers] [-a array] [-d delim] [->
 cd [-L|[-P [-e]] [-@]] [dir]            readarray [-d delim] [-n count] [-O >
 command [-pVv] command [arg ...]        readonly [-aAf] [name[=value] ...] o>
 compgen [-abcdefgjksuv] [-o option] [>  return [n]
 complete [-abcdefgjksuv] [-pr] [-DEI]>  select NAME [in WORDS ... ;] do COMM>
 compopt [-o|+o option] [-DEI] [name .>  set [-abefhkmnptuvxBCHP] [-o option->
 continue [n]                            shift [n]
 coproc [NAME] command [redirections]    shopt [-pqsu] [-o] [optname ...]
 declare [-aAfFgiIlnrtux] [-p] [name[=>  source filename [arguments]
 dirs [-clpv] [+N] [-N]                  suspend [-f]
 disown [-h] [-ar] [jobspec ... | pid >  test [expr]
 echo [-neE] [arg ...]                   time [-p] pipeline
 enable [-a] [-dnps] [-f filename] [na>  times
 eval [arg ...]                          trap [-lp] [[arg] signal_spec ...]
 exec [-cl] [-a name] [command [argume>  true
 exit [n]                                type [-afptP] name [name ...]
 export [-fn] [name[=value] ...] or ex>  typeset [-aAfFgiIlnrtux] [-p] name[=>
 false                                   ulimit [-SHabcdefiklmnpqrstuvxPT] [l>
 fc [-e ename] [-lnr] [first] [last] o>  umask [-p] [-S] [mode]
 fg [job_spec]                           unalias [-a] name [name ...]
 for NAME [in WORDS ... ] ; do COMMAND>  unset [-f] [-v] [-n] [name ...]
 for (( exp1; exp2; exp3 )); do COMMAN>  until COMMANDS; do COMMANDS; done
 function name { COMMANDS ; } or name >  variables - Names and meanings of so>
 getopts optstring name [arg ...]        wait [-fn] [-p var] [id ...]
 hash [-lr] [-p pathname] [-dt] [name >  while COMMANDS; do COMMANDS; done
 help [-dms] [pattern ...]               { COMMANDS ; } parameters.

If you'd like to see a mock of your typer commands as they'll appear in the  document, you can run:  in a properly configured dev environment. Note that this file should not be included in your PRs. The  should be only updated through the git-action, .


