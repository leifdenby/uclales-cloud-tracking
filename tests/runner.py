import os
import subprocess
import shutil

import pytest
import xarray as xr


def run_tracking(
    data_path, base_name, tn_start, tn_end, tracking_type, U_offset=None, move_to=None
):
    if "CMAKE_BINARY_DIR" in os.environ:
        bin_path = os.environ["CMAKE_BINARY_DIR"]
    else:
        print("cmake didn't supply CMAKE_BINARY_DIR")
        bin_path = os.environ["OLDPWD"]

    main_path = os.path.join(bin_path, "main")

    old_cwd = os.getcwd()
    try:
        os.chdir(data_path)

        arg_str = "{} {} {} {} {}".format(
            main_path, base_name, tn_start, tn_end, tracking_type
        )
        if U_offset is not None:
            arg_str += " {} {}".format(*U_offset)

        p = subprocess.Popen(
            arg_str.split(" "), stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )

        stdout, stderr = p.communicate()

        if p.returncode != 0:
            raise Exception(stdout, stderr)

        print(stdout)

        fn_out = os.path.join(data_path, "{}.out.xy.track.nc".format(base_name))

        if move_to is not None:
            if not move_to.endswith(".nc"):
                move_to += ".nc"
            shutil.move(fn_out, move_to)
            fn_out = move_to
        return xr.open_dataset(fn_out, decode_times=False)
    except:
        raise
    finally:
        os.chdir(old_cwd)
