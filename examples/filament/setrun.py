"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np

#------------------------------
def setrun(claw_pkg='amrclaw'):
#------------------------------

    from clawpack.clawutil import data
    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # ---------------
    # Spatial domain:
    # ---------------

    clawdata.num_dim = num_dim

    clawdata.lower[0] = 0          # xlower
    clawdata.upper[0] = 2          # xupper
    clawdata.lower[1] = 0          # ylower
    clawdata.upper[1] = 2          # yupper

    clawdata.num_cells[0] = 64      # mx
    clawdata.num_cells[1] = 64      # my


    clawdata.dt_variable = False
    clawdata.dt_initial = 4e-3

    clawdata.output_style = 3
    # For single grid
    clawdata.output_step_interval = 250
    clawdata.total_steps = 2000

    # clawdata.total_steps = 32000
    # clawdata.output_step_interval = 4000
    clawdata.verbosity = 0


    # ---------------------------
    # AMR parameters (more below)
    # ---------------------------
    amrdata = rundata.amrdata

    amrdata.amr_levels_max = 3

    amrdata.refinement_ratios_x = [4, 4, 2, 2, 2, 2]
    amrdata.refinement_ratios_y = [4, 4, 2, 2, 2, 2]
    amrdata.refinement_ratios_t = [4, 4, 2, 2, 2, 2]

    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.05  # tolerance used in this routine
    amrdata.regrid_interval = 1
    amrdata.regrid_buffer_width  = 2
    amrdata.clustering_cutoff = 0.700000

    # ---------------
    # Size of system:
    # ---------------

    clawdata.num_eqn = 1
    clawdata.num_aux = 3
    clawdata.capa_index = 0


    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.000000


    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data


    # -------------
    # Output times:
    #--------------

#     clawdata.output_style = 3
#
#     if clawdata.output_style==1:
#         # Output ntimes frames at equally spaced times up to tfinal:
#         # Can specify num_output_times = 0 for no output
#         clawdata.num_output_times = 64
#         clawdata.tfinal = 8.0
#         clawdata.output_t0 = True  # output at initial (or restart) time?
#
#     elif clawdata.output_style == 2:
#         # Specify a list or numpy array of output times:
#         # Include t0 if you want output at the initial time.
#         clawdata.output_times =  [0., 0.5, 1.0]
#
#     elif clawdata.output_style == 3:
#         clawdata.output_style = 3
#         # Output every step_interval timesteps over total_steps timesteps:
#         clawdata.output_step_interval = 16000
#         clawdata.total_steps = 2000
#         # clawdata.output_step_interval = 250
#         # clawdata.total_steps = 2000

    clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format = 'ascii'       # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'    # only 'all'
    clawdata.output_aux_components = 'none'  # 'all' or 'none'
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0?


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # --------------
    # Time stepping:
    # --------------

    clawdata.dt_max = 1.000000e+99

    clawdata.cfl_desired = 0.900000
    clawdata.cfl_max = 1.000000
    clawdata.steps_max = 100000


    # ------------------
    # Method to be used:
    # ------------------

    clawdata.order = 2
    clawdata.dimensional_split = 'unsplit'
    clawdata.transverse_waves = 2


    clawdata.num_waves = 1
    clawdata.limiter = ['mc']
    clawdata.use_fwaves = False    # True ==> use f-wave version of algorithms
    clawdata.source_split = 0


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper

    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    rundata.gaugedata.gauges.append([1, 0.4, 0.3, 0., 10.])
    rundata.gaugedata.gauges.append([2, 0.6, 0.3, 0., 10.])


    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------

#     # These are set above.
#     amrdata = rundata.amrdata
#
#     amrdata.amr_levels_max = 1
#
#     amrdata.refinement_ratios_x = [4, 4, 2, 2, 2, 2]
#     amrdata.refinement_ratios_y = [4, 4, 2, 2, 2, 2]
#     amrdata.refinement_ratios_t = [4, 4, 2, 2, 2, 2]

    amrdata.aux_type = ['xleft', 'yleft', 'center']


    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.000000e+00  # Richardson tolerance

    amrdata.verbosity_regrid = 1


    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]


    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting


    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
