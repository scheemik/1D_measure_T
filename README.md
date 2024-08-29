# 1D Measure _T_

by Mikhail Schee

This is a repo to run numerical experiments to calculate the transmission coefficient _T_ for internal waves incident on a vertical stratification structure. Each run can have a different number of mixed layers in the stratification profile and a different value of _kL_, the ratio between the wavelength of the internal wave to the mixed layer depth, _L_. This work was presented in chapter 4 of my thesis, "Thermohaline staircases in the Arctic Ocean: Detection, evolution, and interaction"

This repo runs the experiments in Dedalus 2, a flexible framework for solving PDEs, see [Dedalus Project](https://dedalus-project.readthedocs.io/en/latest/). Because my intention was to run many simulations with different parameters to compare their values of _T_, I have orgainzed the repo as follows:

* Top directory
    * Here, I have the master version of the files I use to run the experiments. The most important of these are `exp_run.sh`, which is what you would launch to run an experiment locally, and `HPC_job_submit.sh`, which is what you would run to submit a job to the high-performance cluster. Each has a number of flags you can envoke to modify which parts of the code will be run. Before launching these, you should make sure you have the parameters set the way you would like them to be in the code files.
* In the `_code_files` directory
    * These are the master copies of the code files. When running in parallel on HPC, directory locking prevents subprocesses from using modifying anything in a directory where a different subprocess is running code. In order to allow many processes to run in parallel, each with a different version of the experiment, I found I needed to make a copy of all the code files for each simulation. The code files in this directory are the ones copies are made from.
    * If you want to change how all future experiments will run, edit the code files here.
* In the `_experiments` directory
    * Each sub-folder under this one represents a different experiment. Each experiment can have any number of specified runs. Each experiment has a copy of the `_code_files` directory, a directory of `_simulations`, and the overall output files.
    * The `_experiments/test_exp/_code_files` directory
        * This is a copy of the `_code_files`. Each simulation will copy these code files to do its run. If you want to change how all simulations for just this experiment (`test_exp`) works, modify these files.
    * `_simulations`
        * Inside this directory, there will be a directory for each simulation run, numbered in order. Within each, there will be a copy of all the code files, snapshots of the simulation data, and the outputs that were specified to make when running the experiment
    * Outputs
        * There will also be outputs here which combine all simulations. Usually, this is a plot of _T_ across all simulations. Sometimes, in order to make these plots, you will need to navigate to the experiment folder and run the code `plot_exp_data.py`. If combining multiple experiments, you will need to copy in the `exp_data.csv` files from each experiment in to one experiment folder and run a modified version of `plot_exp_data.py`. An example can be seen in the `new_128s_many_l` experiment folder.