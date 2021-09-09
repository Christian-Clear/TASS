# TAME

Term Analysis Made Easy (TAME) is designed to make the process of term analysis more user-friendly and promote the standardised storage and transfer of project files.
        
TAME uses a new implementation of the STRANS code to find matching lines in a linelist from a list of user-supplied energy levels. These matched lines are passed to the least-squares optimisation program (LOPT) to calculate optimised energy levels.
        
To download the latest TAME version or submit bug reports and improvement suggestions, please visit the project homepage. 

## Installation

tame is written entirely in Python and runs through the Python interpreter. As python is not a compiled language, libraries that \tame depends on must be installed on each computer that runs it. Luckily, there is an easy way to accomplish this! 

The popular scientific toolkit for python, Anaconda, is able to run programs within specific environments. Included in the tame files is the YAML file, tame.yml, which contains all of the information Anaconda (specifically the package manager Conda) needs to create an environment with exactly the correct packages (and versions of packages) for tame to run correctly.

To create the Anaconda environment, first download and install Anaconda. Then navigate to the tame files and create the tame environment:

`(base) user@host: cd <main tame directory /resources>`

`(base) user@host: conda create env --file tame.yml`

this will create the TASS environment in your main Anaconda /envs folder. If you wish the tame environment to be saved elsewhere, then you should use --prefix tag:

`(base) user@host: conda create env --file tame.yml --prefix ./<directory>`

The creation of the environment may take some time (there are a lot of libraries to be installed!).

## Running TAME

Once the tame environment has been created, you need to activate it:

`(base) user@host: conda activate tame`

`(tame) user@host: `

(note the change from base to tame in the brackets of the terminal).

You can now run tame from within the tame environment (you will need to navigate to the main tame folder):

`(tame) user@host: cd <main tame directory>`

`(tame) user@host: python tame.py`
