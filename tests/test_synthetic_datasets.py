from utils import BackgroundState, FakeCloudState, FakeSimulation

from runner import run_tracking

import os

def test_single_consecutive_clouds():
    background = BackgroundState()

    clouds_list = []
    t0 = 300.
    n_clouds = 2
    for n in range(n_clouds):
        cloud = FakeCloudState(background=background, t0=t0, v=(1., 2.))
        t0 += cloud.t_max + 120.
        clouds_list.append(cloud)

    fake_sim = FakeSimulation(clouds_list=clouds_list, background=background)

    with fake_sim.output(keep_files=True) as (path, ds_input):
        print(path)
        ds_track = run_tracking(data_path=path, base_name='testdata',
                                tn_start=1, tn_end=len(ds_input.time),
                                tracking_type='cloud,core')

        assert ds_track.smcloudid.max() == n_clouds
        ds_track.close()


        # ensure we can track not starting from the first timestep
        ds_track = run_tracking(data_path=path, base_name='testdata',
                                tn_start=2, tn_end=len(ds_input.time),
                                tracking_type='cloud,core')

        assert ds_track.smcloudid.max() == n_clouds
        ds_track.close()


# def test_one(ds_track):

def foobar():
    v0 = (1., 1.)

    clouds_list = [
       FakeCloudState(background=background, t0=400, v=v0),
       FakeCloudState(background=background, p0=(500., 0.), v=v0),
       FakeCloudState(background=background, p0=(-500., 500.), v=v0), 
       FakeCloudState(background=background, t0=400., p0=(1000., 1000.), v=v0),    
       FakeCloudState(background=background, t0=800., p0=(1000., 1000.), v=v0),    

    ]

    fake_sim = FakeSimulation(clouds_list=clouds_list, background=background)

    fake_sim(t=2400.).lwp.plot()
#clouds(800.).z_top.plot()
#clouds_list[0](x_, y_, 100).core.plot()


# In[431]:




# In[432]:


