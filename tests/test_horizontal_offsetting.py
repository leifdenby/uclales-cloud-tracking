import pytest
import numpy as np

from utils import BackgroundState, FakeCloudState, FakeSimulation
from runner import run_tracking


def _build_offsetting_test_sim(is_offset, offset_dim=None):
    r_cloud = 300.0
    dt = 120

    if is_offset:
        if offset_dim == "x":
            U_offset = (2 * r_cloud / dt, 0.0)
        elif offset_dim == "y":
            U_offset = (0.0, 2 * r_cloud / dt)
        else:
            raise NotImplementedError(offset_dim)
    else:
        U_offset = (
            0.0,
            0.0,
        )

    background = BackgroundState()

    clouds_list = [
        FakeCloudState(
            background=background,
            p0=(500.0, 0.0),
            v=U_offset,
            r_max=r_cloud,
            r_min=r_cloud,
        ),
    ]

    return (
        FakeSimulation(clouds_list=clouds_list, background=background),
        dt,
        U_offset,
    )


fast_moving_cloud_sims = [
    _build_offsetting_test_sim(is_offset=True, offset_dim="x"),
    _build_offsetting_test_sim(is_offset=True, offset_dim="y"),
]


@pytest.mark.parametrize("fake_sim,dt,U_offset", fast_moving_cloud_sims)
def test_single_fast_moving_cloud_number(fake_sim, dt, U_offset):
    """
    To test correction for Galilean transform (or wind in simulation) make some
    clouds that will move far enough to not have spatial overlap between
    timesteps and check that there is still only one cloud when applying the
    offset
    """

    bn = "test_moving_clouds_{}_{}".format(*U_offset)

    with fake_sim.output(dt_sim=dt, base_name=bn) as (path, ds_input):
        # first without the offset applied during tracking
        ds_track_no_offset = run_tracking(
            data_path=path,
            base_name=ds_input.base_name,
            tn_start=1,
            tn_end=len(ds_input.time),
            tracking_type="cloud,core",
            U_offset=None,
            move_to="tracking_no_offset",
        )

        # with no offset there will be more than one cloud
        assert ds_track_no_offset.smcloudid.max() != 1

        # now with the offset applied
        ds_track_offset = run_tracking(
            data_path=path,
            base_name=ds_input.base_name,
            tn_start=1,
            tn_end=len(ds_input.time),
            tracking_type="cloud,core",
            U_offset=U_offset,
            move_to="tracking_offset",
        )

        # with the offset applied there should only be one cloud
        assert ds_track_offset.smcloudid.max() == 1

        # although the cloud mask should now only include one cloud it should
        # still be in the same xy-locations as before
        m_iscloud_1 = ds_track_no_offset.nrcloud.fillna(0) > 0
        m_iscloud_2 = ds_track_offset.nrcloud.fillna(0) == 1
        assert np.all(m_iscloud_1 == m_iscloud_2)

        ds_track_no_offset.close()
        ds_track_offset.close()
