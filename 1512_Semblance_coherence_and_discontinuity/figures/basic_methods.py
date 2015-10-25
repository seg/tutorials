import numpy as np
import scipy.ndimage as ndimage
import scipy.signal

def bahorich_coherence(data, zwin):
    ni, nj, nk = data.shape
    out = np.zeros_like(data)
    padded = np.pad(data, ((0, 0), (0, 0), (zwin//2, zwin//2)), mode='reflect')

    for i, j, k in np.ndindex(ni - 1, nj - 1, nk - 1):
        center_trace = data[i,j,:]
        center_std = center_trace.std()
        x_trace = padded[i+1, j, k:k+zwin]
        y_trace = padded[i, j+1, k:k+zwin]

        xcor = np.correlate(center_trace, x_trace)
        ycor = np.correlate(center_trace, y_trace)

        px = xcor.max() / (xcor.size * center_std * x_trace.std())
        py = ycor.max() / (ycor.size * center_std * y_trace.std())
        out[i,j,k] = np.sqrt(px * py)

    return out

def moving_window(data, window, func):
    wrapped = lambda region: func(region.reshape(window))
    return ndimage.generic_filter(data, wrapped, window, mode='reflect')

def marfurt_semblance(region):
    region = region.reshape(-1, region.shape[-1])
    ntraces, nsamples = region.shape

    cov = region.dot(region.T)
    return cov.sum() / cov.diagonal().sum() / ntraces

def semblance2(region):
    region = region.reshape(-1, region.shape[-1])
    ntraces, nsamples = region.shape

    square_of_sums = np.sum(region, axis=0)**2
    sum_of_squares = np.sum(region**2, axis=0)
    return square_of_sums.sum() / sum_of_squares.sum() / ntraces

def eig(region):
    region = region.reshape(-1, region.shape[-1])

    cov = region.dot(region.T)
    vals = np.linalg.eigvalsh(cov)
    return vals.max() / vals.sum()

def complex_semblance(region):
    region = region.reshape(-1, region.shape[-1])
    ntraces, nsamples = region.shape

    region = scipy.signal.hilbert(region)
    region = np.hstack([region.real, region.imag])
    cov = region.dot(region.T)
    return np.abs(cov.sum() / cov.diagonal().sum()) / ntraces

def complex_eig(region):
    region = region.reshape(-1, region.shape[-1])

    region = scipy.signal.hilbert(region)
    region = np.hstack([region.real, region.imag])

    cov = region.dot(region.T)
    vals = np.linalg.eigvals(cov)
    return np.abs(vals.max() / vals.sum())

def flatten(data, surface, window):
    surface = ndimage.gaussian_filter(surface.astype(float), 3)

    ni, nj, nk = data.shape
    ik = np.arange(nk)
    out_ik = np.arange(window) - window // 2

    out = np.zeros((ni, nj, window))
    for i, j in np.ndindex(ni, nj):
        trace = data[i,j,:]
        k = surface[i, j]
        shifted = np.interp(out_ik + k, ik, trace)

        out[i,j,:] = shifted

    return out

def unflatten(data, surface, orig_shape):
    out = np.zeros(orig_shape)
    surface = np.clip(surface, 0, orig_shape[-1] - 1)

    win = data.shape[-1] // 2
    for i, j in np.ndindex(orig_shape[0], orig_shape[1]):
        k = surface[i,j]

        outmin, outmax = max(0, k - win), min(orig_shape[-1], k + win + 1)
        inmin, inmax = outmin - (k - win), k + win + 1 - outmax
        inmax = data.shape[-1] - abs(inmax)

        out[i, j, outmin:outmax] = data[i, j, inmin:inmax]

    return out
