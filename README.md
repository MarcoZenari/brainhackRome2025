## Download the Data from OpenNeuro
Change working directory to the repository folder and run the following command
```bash
aws s3 sync --no-sign-request s3://openneuro.org/ds005237 ./data \
  --exclude "*" \
  --include "fMRI_timeseries_clean_denoised_GSR_parcellated/*/*" \
  --include "motion_FD/*" \
  --include "phenotype/*"
```
then we are going to move some files to make the code run easier
```bash
mv ./data/phenotype/demos.tsv ./data/
mv ./data/phenotype/notes.tsv ./data/
```
```
aws s3 sync --no-sign-request s3://openneuro.org/ds005237 ./data/events \
  --exclude "*" \
  --include "sub-*/func/*events.tsv"
```
## Install ssm
Activate the environment you want to use, then:
```bash
git submodule add https://github.com/lindermanlab/ssm.git submodules/ssm
cd submodules/ssm
pip install numpy cython
pip install -e .
```
