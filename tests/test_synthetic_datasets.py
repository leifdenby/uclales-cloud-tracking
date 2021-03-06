import pytest

from utils import BackgroundState, FakeCloudState, FakeSimulation
from runner import run_tracking


def test_single_consecutive_clouds():
    background = BackgroundState()

    clouds_list = []
    t0 = 300.0
    n_clouds = 2
    for n in range(n_clouds):
        cloud = FakeCloudState(background=background, t0=t0, v=(1.0, 2.0))
        t0 += cloud.t_max + 120.0
        clouds_list.append(cloud)

    fake_sim = FakeSimulation(clouds_list=clouds_list, background=background)

    with fake_sim.output() as (path, ds_input):
        ds_track = run_tracking(
            data_path=path,
            base_name="testdata",
            tn_start=1,
            tn_end=len(ds_input.time),
            tracking_type="cloud,core",
        )

        assert ds_track.smcloudid.max() == n_clouds
        ds_track.close()

        # test bug hasn't been fixed
        # (https://github.com/leifdenby/uclales-cloud-tracking/issues/4)
        with pytest.raises(Exception):
            # ensure we can track not starting from the first timestep
            ds_track = run_tracking(
                data_path=path,
                base_name="testdata",
                tn_start=2,
                tn_end=len(ds_input.time),
                tracking_type="cloud,core",
            )

            assert ds_track.smcloudid.max() == n_clouds
            ds_track.close()


def test_two_merging_clouds():
    """
    The tracking code merges two clouds into one if their cores overlap at some
    point in their lifetime
    """
    background = BackgroundState()

    clouds_list = [
        FakeCloudState(background=background, v=(1.0, 1.0)),
        FakeCloudState(background=background, p0=(1000.0, 1000.0)),
    ]

    fake_sim = FakeSimulation(clouds_list=clouds_list, background=background)

    with fake_sim.output() as (path, ds_input):
        ds_track = run_tracking(
            data_path=path,
            base_name="testdata",
            tn_start=1,
            tn_end=len(ds_input.time),
            tracking_type="cloud,core",
        )

        assert ds_track.smcloudid.max() == 1
        ds_track.close()


def test_two_splitting_clouds():
    """
    The tracking code will treat to clouds that split apart as one if their
    cores overlap at some point in their lifetime
    """
    background = BackgroundState()

    clouds_list = [
        FakeCloudState(background=background, r_min=300, v=(1.0, 1.0)),
        FakeCloudState(background=background, r_min=300, p0=(1000.0, 1000.0)),
    ]

    fake_sim = FakeSimulation(clouds_list=clouds_list, background=background)

    with fake_sim.output() as (path, ds_input):
        ds_track = run_tracking(
            data_path=path,
            base_name="testdata",
            tn_start=1,
            tn_end=len(ds_input.time),
            tracking_type="cloud,core",
        )

        assert ds_track.smcloudid.max() == 1
        ds_track.close()
