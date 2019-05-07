from pixell.lensing import lens_map_curved
from pixell import curvedsky, enmap

import numpy as np

def rand_map(shape, wcs, ps_lensinput, lmax=None, maplmax=None, dtype=np.float64, seed=None, phi_seed=None, oversample=2.0, spin=[0,2], output="l", geodesic=True, verbose=False, delta_theta=None):
	from pixell import curvedsky, sharp
	ctype   = np.result_type(dtype,0j)
	# Restrict to target number of components
	oshape  = shape[-3:]
	if len(oshape) == 2: shape = (1,)+tuple(shape)
	ncomp   = shape[-3]
	ps_lensinput = ps_lensinput[:1+ncomp,:1+ncomp]
	# First draw a random lensing field, and use it to compute the undeflected positions
	if verbose: print("Generating alms")
	if phi_seed is None:
		alm, ainfo = curvedsky.rand_alm(ps_lensinput, lmax=lmax, seed=seed, dtype=ctype, return_ainfo=True)
	else:
		# We want separate seeds for cmb and phi. This means we have to do things a bit more manually
		wps, ainfo = curvedsky.prepare_ps(ps_lensinput, lmax=lmax)
		alm = np.empty([1+ncomp,ainfo.nelem],ctype)
		curvedsky.rand_alm_white(ainfo, alm=alm[:1], seed=phi_seed)
		curvedsky.rand_alm_white(ainfo, alm=alm[1:], seed=seed)
		ps12 = enmap.multi_pow(wps, 0.5)
		ainfo.lmul(alm, (ps12/2**0.5).astype(dtype), alm)
		alm[:,:ainfo.lmax].imag  = 0
		alm[:,:ainfo.lmax].real *= 2**0.5
		del wps, ps12
		
	phi_alm, cmb_alm = alm[0], alm[1:]
	del alm
	# Truncate alm if we want a smoother map. In taylens, it was necessary to truncate
	# to a lower lmax for the map than for phi, to avoid aliasing. The appropriate lmax
	# for the cmb was the one that fits the resolution. FIXME: Can't slice alm this way.
	#if maplmax: cmb_alm = cmb_alm[:,:maplmax]
	return lens_map_curved(shape=shape, wcs=wcs, phi_alm=phi_alm,
						   cmb_alm=cmb_alm, phi_ainfo=ainfo, maplmax=maplmax,
						   dtype=dtype, oversample=oversample, spin=spin,
						   output=output, geodesic=geodesic, verbose=verbose,
						   delta_theta=delta_theta)
