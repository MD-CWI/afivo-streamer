# Parameter study

## Usage

1. edit the file `generate_parameters.py`, in particular the number of samples, the "base" .cfg file to use, and the parameters that should be varied
2. run `generate_parameters.py` which will produce `commands.txt`
3. run all simulations by executing `run_study.sh`, which requires GNU parallel to be installed. You can modify `run_study.sh` beforehand to change the number of simultaneous runs and threads
4. execute `store_npz_data.sh` (after potentially modifying it) to convert .silo files to .npz files for each variable of interest
5. run `create_dataset.sh`, which will execute `create_dataset.py` multiple times to create a dataset containing different variables

Optionally, you can look at some examples in the dataset with `plot_dataset.py`


