# Two-phase VoF with a static IBM cylinder

This experiment compares two otherwise identical Stokes II wave cases:

* `stokesII-control`: IBM disabled;
* `stokesII-static-ibm`: a fixed SDF cylinder with four MDF passes.

The wave period and ramp time are both `3 s`.  Each case runs to `9 s`,
covering the ramp and two additional wave periods.  The maximum time
step is `0.01 s`.

Run-time cases and output are placed under:

`~/fastscratch/ichiBoundaryFoam-2phase-static-261104/`

This is a first coupled verification.  It checks phase boundedness,
Courant numbers, volume drift, IBM slip, projection divergence, and
changes relative to the no-IBM wave control.

