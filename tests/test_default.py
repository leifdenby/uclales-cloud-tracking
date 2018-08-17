import os
import subprocess

import pytest
import xarray as xr


@pytest.fixture
def ds_track(request):
    bin_path = os.environ['OLDPWD']
    main_path = os.path.join(bin_path, 'main')
    test_path = os.path.dirname(os.path.abspath(__file__))
    testdata_path = os.path.join(test_path, "testdata")

    print(">>", os.environ['OLDPWD'], os.environ['PWD'])

    os.chdir(testdata_path)

    base_name = "testdata"
    arg_str = "{} {} 1 5 {}".format(
        main_path, base_name, request.param
    )
    p = subprocess.Popen(arg_str.split(" "), stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)

    stdout, stderr = p.communicate()

    if p.returncode != 0:
        raise Exception(stdout, stderr)

    fn_out = os.path.join(testdata_path, '{}.out.xy.track.nc'.format(base_name))

    return xr.open_dataset(fn_out, decode_times=False)

@pytest.mark.parametrize("ds_track", ["cloud", ], indirect=["ds_track",])
def test_one(ds_track):
    assert ds_track.smcloudid.max() == 6
