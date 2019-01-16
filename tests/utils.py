# coding: utf-8
"""
Utility to generate synthetic dataset for testing with parameterised fields of
clouds
"""

import xarray as xr
import numpy as np

import os
import shutil


lx, ly = 4000., 4000.,
z_min, z_max = 600., 2000.,
dx = dy = 25.


DEFAULT_STATE_VARS = "z_top z_base core lwp".split(" ")
UNITS = dict(z_top='m', z_base='m', core='K', lwp='m g/kg')

class StateBaseClass(object):   
    def __call__(self, t):
        da_ = []
        
        for var_name in self.state_vars:
            fn = getattr(self, '_get_{}'.format(var_name))
            values = fn(t)
            da = xr.DataArray(
                values,
                coords=dict(x=self.x, y=self.y), dims=('y', 'x'), 
                attrs=dict(time=t, units=UNITS[var_name]),
                name=var_name
            )
            da_.append(da)
            
        ds = xr.merge(da_)
        return ds.expand_dims('time')    

class BackgroundState(StateBaseClass):
    def __init__(self, lx=lx, ly=ly, state_vars=DEFAULT_STATE_VARS):
        self.state_vars = state_vars
        self.x = np.linspace(-lx/2., lx/2, int(lx/dx)+1)
        self.y = np.linspace(-ly/2, ly/2, int(ly/dy)+1)
        self.nx = len(self.x)
        self.ny = len(self.y)

    def _get_z_top(self, t):
        return np.zeros((self.nx, self.ny))

    def _get_z_base(self, t):
        return np.zeros((self.nx, self.ny))

    def _get_core(self, t):
        return np.zeros((self.nx, self.ny))

    def _get_lwp(self, t):
        return np.zeros((self.nx, self.ny))

class FakeCloudState(StateBaseClass):
    def __init__(self, background, p0=(0., 0.), v=(0., 0.), t0=0,
                 tau=60.*15., core_max=2.0, lwp_max=10., r_min=100., 
                 r_max=500.):
        """
        Represents a single fake cloud

        p0: (x, y) origin [m]
        v: (u, v) velocity of cloud [m/s]
        t0: birth time of cloud [s]
        tau: growth timescale [s]
        core_max: max core potential temperature anamoly [K]
        lwp_max: max integrated liquid water path [m g/kg]
        """
        self.p0 = p0
        self.t0 = t0
        self.v = v
        self.x = background.x
        self.y = background.y
        self.state_vars = background.state_vars
        self.nx = len(self.x)
        self.ny = len(self.y)

        self.tau = tau
        self.lwp_max = lwp_max
        self.core_max = core_max

        self.t_max = self.t0 + 2*self.tau
        self.r_min = r_min
        self.r_max = r_max

    def r(self, t):
        @np.vectorize
        def r_fn(t):
            t_ = t - self.t0
            if t_ < 0.:
                return 0.0
            elif t < self.t_max:
                dr = self.r_max - self.r_min
                return self.r_min + dr*np.min([1., t_/self.tau])
            else:
                return 0.0
        
        return r_fn(t)
    
    @staticmethod
    def _wrap_range(l, p, p_min):
        return p - l*np.round((p-l/2.-p_min)/l)
    
    def _get_offset_grid(self, t):
        xx, yy = np.meshgrid(self.x, self.y,)
        
        x_c = self.p0[0] - self.v[0]*(t-self.t0)    
        y_c = self.p0[1] - self.v[1]*(t-self.t0)
        
        xx -= x_c
        yy -= y_c
        
        xx = self._wrap_range(lx, xx, self.x.min())
        yy = self._wrap_range(ly, yy, self.y.min())

        return xx, yy

    
    def _get_mask(self, t):
        r = self.r(t)

        xx, yy = self._get_offset_grid(t)
        
        dist = np.sqrt(xx**2. + yy**2.)*np.sign(t - self.t0)
        m = np.logical_and(dist < r, r > 0.0)
        return m, r, dist
    
    def _get_z_base(self, t):
        m, _, _ = self._get_mask(t)
        return xr.where(m, z_min, np.nan)
    
    def _get_z_top(self, t):
        @np.vectorize
        def z_top_fn(t):
            t_ = t - self.t0
            if t_ < 0.:
                return np.nan
            else:
                return z_min + (z_max - z_min)*np.max([0., np.min([t_/self.tau, 1.0])])
            
        z_top = z_top_fn(t)
        m, _, _ = self._get_mask(t)
    
        z_ = np.ones((self.nx, self.ny))*np.nan
        z_[m] = z_top
        
        return z_

    def _get_core(self, t):
        m, r, dist = self._get_mask(t)
        return xr.where(m, self.core_max*(1. - dist/r), np.nan)
    
    def _get_lwp(self, t):
        m, r, dist = self._get_mask(t)
        return xr.where(np.logical_and(dist < r, r > 0.0), self.lwp_max*(1. - dist/r), np.nan)

class FakeSimulation:
    def __init__(self, clouds_list, background):
        self.clouds_list = clouds_list
        self.background = background
        
    @staticmethod
    def _combine_clouds(clds):
        for n, c in enumerate(clds):
            c.expand_dims('cloud_id')
            c.coords['cloud_id'] = n+1

        ds_clouds = xr.concat(clds, dim='cloud_id')

        def reduce_da(da):
            if da.name == 'z_top':
                return da.max(dim='cloud_id', keep_attrs=True)
            elif da.name == 'z_base':
                return da.min(dim='cloud_id', keep_attrs=True)
            else:
                return da.sum(dim='cloud_id', keep_attrs=True)

        return ds_clouds.apply(reduce_da)
        
    def __call__(self, t):
        ds_clouds = self._combine_clouds([c(t) for c in self.clouds_list])
        ds_background = self.background(t)

        def merge_cloud_and_background_fields(da):
            return da.where(~da.isnull(), ds_background[da.name])

        return ds_clouds.apply(merge_cloud_and_background_fields)

    @property
    def t_max(self):
        return np.max([c.t_max for c in self.clouds_list])

    def generate_output(self, dt=2*60):
        t_ = np.arange(0., self.t_max + 120., dt)
        ds = xr.concat([self(t) for t in t_], dim='time')
        ds.coords['time'] = t_
        ds.coords['time'].attrs['units'] = 'seconds since 1970-01-01 00:00:00'

        return ds.rename(dict(x='xt', y='yt', z_top='cldtop', z_base='cldbase'))

    def output(self, dt_sim=2*60, keep_files=False):
        return SyntheticDatasetGenerator(sim=self, dt_sim=dt_sim,
                                         keep_files=keep_files)


class SyntheticDatasetGenerator:
    BASE_NAME = 'testdata'
    def __init__(self, sim, dt_sim, keep_files=False):
        self.sim = sim
        self.dt_sim = dt_sim
        self.keep_files = keep_files

    def __enter__(self):
        ds_input = self.sim.generate_output(dt=self.dt_sim)

        path = os.path.join(os.getcwd(), self.BASE_NAME)
        if not os.path.exists(path):
            os.makedirs(path)
        fn_base = "{}.out.xy.{}.nc"
        for v in ds_input.data_vars:
            fn = os.path.join(path, fn_base.format(self.BASE_NAME, v))
            ds_input[v].to_netcdf(fn)

        self.path = path
        return path, ds_input

    def __exit__(self, *args, **kwargs):
        if not self.keep_files:
            shutil.rmtree(self.path)


if __name__ == "__main__":
    background = BackgroundState()

    clouds_list = []
    t0 = 300.
    n_clouds = 2
    for n in range(n_clouds):
        cloud = FakeCloudState(background=background, t0=t0, v=(1., 2.))
        t0 += cloud.t_max + 120.
        clouds_list.append(cloud)

    fake_sim = FakeSimulation(clouds_list=clouds_list, background=background)

    path = os.path.join(os.getcwd(), 'testdata')

    if not os.path.exists(path):
        os.makedirs(path)

    fake_sim.generate_output(path=path)
