# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

# Import abipy modules
from abipy import abilab
import abipy.data as abidata
from abipy import flowtk
import os
import sys

def main():

    structure = abidata.structure_from_ucell("GaAs")
    print(structure.abi_string)
    pseudos = abidata.pseudos("31ga.pspnc","33as.pspnc")
    band_flow = build_flow(structure,pseudos)
    print(band_flow[0][2].input)
    band_flow.make_scheduler().start()
    
def make_scf_nscf_inputs(structure,pseudos,options=None):
    multi = abilab.MultiDataset(structure=structure, 
                                pseudos=pseudos, ndtset=4)
    # Global variables
    ecut = 15
    multi.set_vars(
        ecut=ecut,
        nband=20,
        iomode=3,
        timopt=-1,
        diemac=14.0,
        nbdbuf=2
    )

    shiftk = [0,0,0]
    # Dataset 0 (GS run)
    multi[0].set_kmesh(ngkpt=[12, 12, 12], shiftk=shiftk)
    multi[0].set_vars(toldfe=1e-8)


    # Dataset 1 (NSCF run with k-path)
    multi[1].set_vars(iscf=-2,
                      tolwfr=1e-12,
                      nband=20,
                      )
    multi[1].set_kpath(ndivsm=8)
    
    nscf_kmesh = [32,32,32]
    # Dataset 2 (NSCF run)
    multi[2].set_kmesh(ngkpt=nscf_kmesh, shiftk=shiftk)
    multi[2].set_vars(tolwfr=1e-16, iscf=-2)
    
    multi[3].set_kmesh(ngkpt=nscf_kmesh, shiftk=shiftk)
    multi[3].set_vars(
        optdriver=8,
        wfk_task='"wfk_optics"',
        autoparal=0,
        max_ncpus=0,
        iscf=-2
    )
    # Generate two input files for the GS and the NSCF run
    scf_input, bs_input, nscf_input, optics_input = multi.split_datasets()
    return scf_input, bs_input, nscf_input, optics_input

def build_flow(structure,pseudos):
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Get the SCF and the NSCF input.
    scf_input, bs_input, nscf_input, optics_input = make_scf_nscf_inputs(structure,pseudos)

    # Build the flow.
    flow = flowtk.Flow(workdir=workdir)
    flow.register_task(scf_input)
    flow.register_task(bs_input,deps={flow[0][0]: "DEN"},append=True)
    flow.register_task(nscf_input,deps={flow[0][0]: "DEN"},append=True)
    flow.register_task(optics_input,deps={flow[0][2]: "WFK"},append=True)
    return flow
    #return flowtk.bandstructure_flow(workdir, scf_input, nscf_input)

if __name__ == "__main__":
    main()
