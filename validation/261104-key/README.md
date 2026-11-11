# Four key IBM validations

This suite is the first, short screening stage for:

1. spatial and temporal convergence;
2. comparison with the body-fitted cylinder;
3. preliminary IBM load reconstruction; and
4. conservation during prescribed body motion.

The generated cases and all run-time output live below
`~/fastscratch/ichiBoundaryFoam-validation-261104/`.  This directory
contains only the reproducible setup, scheduler, and analysis tools.

The screening runs end at `t = 0.5 s`, except for the translating case,
which ends at `t = 5 s` so that the body moves through one medium-grid
cell.  Passing this stage does not establish a converged physical wake.
The next stage must run selected cases for several cylinder convective
times, `D/U`, before comparing mean drag, lift, and shedding frequency.

Run:

```bash
./prepare-suite.sh
cd ~/fastscratch/ichiBoundaryFoam-validation-261104
sbatch sbatch-screening.sh
```

The array contains the following cases:

* `space-coarse`, `space-medium`, `space-fine`: fixed
  `deltaT = 0.00125 s`;
* `time-dt005`, `time-dt0025`, `time-dt00125`: medium mesh with
  fixed time steps;
* `body-fitted-medium`: body-fitted reference;
* `moving-medium`: prescribed translation at `0.002 m/s` for `5 s`.

