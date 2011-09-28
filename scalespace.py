import numpy as np
from numpy import *
import scipy.signal


def companion(a):
    a = np.atleast_1d(a)

    if a.ndim != 1:
        raise ValueError("Incorrect shape for `a`.  `a` must be one-dimensional.")

    if a.size < 2:
        raise ValueError("The length of `a` must be at least 2.")

    if a[0] == 0:
        raise ValueError("The first coefficient in `a` must not be zero.")

    first_row = -a[1:]/(1.0*a[0])
    n = a.size
    c = np.zeros((n-1, n-1), dtype=first_row.dtype)
    c[0] = first_row
    c[range(1,n-1), range(0, n-2)] = 1
    return c


def axis_slice(a, start=None, stop=None, step=None, axis=-1):
    a_slice = [slice(None)] * a.ndim
    a_slice[axis] = slice(start, stop, step)
    b = a[a_slice]
    return b


def axis_reverse(a, axis=-1):
    return axis_slice(a, step=-1, axis=axis)


def odd_ext(x, n, axis=-1):
    if n < 1:
        return x
    if n > x.shape[axis] - 1:
        raise ValueError(("The extension length n (%d) is too big. " +
                         "It must not exceed x.shape[axis]-1, which is %d.")
                         % (n, x.shape[axis] - 1))
    left_end = axis_slice(x, start=0, stop=1, axis=axis)
    left_ext = axis_slice(x, start=n, stop=0, step=-1, axis=axis)
    right_end = axis_slice(x, start=-1, axis=axis)
    right_ext = axis_slice(x, start=-2, stop=-(n + 2), step=-1, axis=axis)
    ext = np.concatenate((2 * left_end - left_ext,
                          x,
                          2 * right_end - right_ext),
                         axis=axis)
    return ext


def even_ext(x, n, axis=-1):
    if n < 1:
        return x
    if n > x.shape[axis] - 1:
        raise ValueError(("The extension length n (%d) is too big. " +
                         "It must not exceed x.shape[axis]-1, which is %d.")
                         % (n, x.shape[axis] - 1))
    left_ext = axis_slice(x, start=n, stop=0, step=-1, axis=axis)
    right_ext = axis_slice(x, start=-2, stop=-(n + 2), step=-1, axis=axis)
    ext = np.concatenate((left_ext,
                          x,
                          right_ext),
                         axis=axis)
    return ext


def const_ext(x, n, axis=-1):
    if n < 1:
        return x
    left_end = axis_slice(x, start=0, stop=1, axis=axis)
    ones_shape = [1] * x.ndim
    ones_shape[axis] = n
    ones = np.ones(ones_shape, dtype=x.dtype)
    left_ext = ones * left_end
    right_end = axis_slice(x, start=-1, axis=axis)
    right_ext = ones * right_end
    ext = np.concatenate((left_ext,
                          x,
                          right_ext),
                         axis=axis)
    return ext



def new_lfilter_zi(b, a):
    # We could use scipy.signal.normalize, but it uses warnings in
    # cases where a ValueError is more appropriate, and it allows
    # b to be 2D.
    b = np.atleast_1d(b)
    if b.ndim != 1:
        raise ValueError("Numerator b must be rank 1.")
    a = np.atleast_1d(a)
    if a.ndim != 1:
        raise ValueError("Denominator a must be rank 1.")

    while len(a) > 1 and a[0] == 0.0:
        a = a[1:]
    if a.size < 1:
        raise ValueError("There must be at least one nonzero `a` coefficient.")

    if a[0] != 1.0:
        # Normalize the coefficients so a[0] == 1.
        a = a / a[0]
        b = b / a[0]

    n = max(len(a), len(b))

    # Pad a or b with zeros so they are the same length.
    if len(a) < n:
        a = np.r_[a, np.zeros(n - len(a))]
    elif len(b) < n:
        b = np.r_[b, np.zeros(n - len(b))]

    IminusA = np.eye(n - 1) - companion(a).T
    B = b[1:] - a[1:] * b[0]
    # Solve zi = A*zi + B
    zi = np.linalg.solve(IminusA, B)

    # For future reference: we could also use the following
    # explicit formulas to solve the linear system:
    #
    # zi = np.zeros(n - 1)
    # zi[0] = B.sum() / IminusA[:,0].sum()
    # asum = 1.0
    # csum = 0.0
    # for k in range(1,n-1):
    #     asum += a[k]
    #     csum += b[k] - a[k]*b[0]
    #     zi[k] = asum*zi[0] - csum

    return zi

def new_filtfilt(b, a, x, axis=-1, padtype='odd', padlen=None):
    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype.  padtype must "
                         "be 'even', 'odd', 'constant', or None.") %
                            padtype)

    b = np.asarray(b)
    a = np.asarray(a)
    x = np.asarray(x)

    ntaps = max(len(a), len(b))

    if padtype is None:
        padlen = 0

    if padlen is None:
        # Original padding; preserved for backwards compatibility.
        edge = ntaps * 3
    else:
        edge = padlen

    # x's 'axis' dimension must be bigger than edge.
    if x.shape[axis] <= edge:
        raise ValueError("The length of the input vector x must be at least "
                         "padlen, which is %d." % edge)

    if padtype is not None and edge > 0:
        # Make an extension of length `edge` at each
        # end of the input array.
        if padtype == 'even':
            ext = even_ext(x, edge, axis=axis)
        elif padtype == 'odd':
            ext = odd_ext(x, edge, axis=axis)
        else:
            ext = const_ext(x, edge, axis=axis)
    else:
        ext = x

    # Get the steady state of the filter's step response.
    zi = new_lfilter_zi(b, a)

    # Reshape zi and create x0 so that zi*x0 broadcasts
    # to the correct value for the 'zi' keyword argument
    # to lfilter.
    zi_shape = [1] * x.ndim
    zi_shape[axis] = zi.size
    zi = np.reshape(zi, zi_shape)
    x0 = axis_slice(ext, stop=1, axis=axis)

    # Forward filter.
    (y, zf) = lfilter(b, a, ext, zi=zi * x0)

    # Backward filter.
    # Create y0 so zi*y0 broadcasts appropriately.
    y0 = axis_slice(y, start=-1, axis=axis)
    (y, zf) = lfilter(b, a, axis_reverse(y, axis=axis), zi=zi * y0)

    # Reverse y.
    y = axis_reverse(y, axis=axis)

    if edge > 0:
        # Slice the actual signal from the extended signal.
        y = axis_slice(y, start=edge, stop=-edge, axis=axis)

    return y





def extend_odd(d, k):
    ld = len(d)
    while len(d) < ld + 2*k:
        d = concatenate([flipud(2*d[0]-d[1:k+1]),
                         d,
                         flipud(2*d[-1]-d[-k-1:-1])])
    if len(d) > ld + 2*k:
        extra = len(d) - (ld + 2*k)
        assert extra % 2 == 0
        d = d[extra/2:-extra/2]
    assert len(d) == ld + 2*k
    return d

def gaussian(data, sigma):
    k = int(ceil(sigma * 2.5))
    d = extend_odd(data, k)
    w = scipy.signal.gaussian(2*k+1, sigma)
    w = w/sum(w)
    print sum(w)
    smoothed = convolve(d, w, 'valid')
    assert len(smoothed) == len(data)
    return smoothed


def gaussian2(data, sigma):
    if sigma >= 2.5:
        q = 0.98711 * sigma - 0.96330
    elif 0.5 <= sigma < 2.5:
        q = 3.97156 - 4.14554 * sqrt(1 - 0.26891 * sigma)
    else:
        assert False

    b0 = 1.57825 + 2.44413*q + 1.4281*q**2 + 0.422205*q**3
    b1 = 2.44413*q + 2.85619*q**2 + 1.26661*q**3
    b2 = -(1.4281*q**2 + 1.26661*q**3)
    b3 = 0.422205*q**3
    B = 1 - ((b1 + b2 + b3)/b0)

    fa = array([b0, -b1, -b2, -b3])
    fb = array([B*b0])
    
    smoothed = new_filtfilt(fb, fa, data)
    assert len(smoothed) == len(data)
    return smoothed


def downsample(d):
    return d[::2]

def upsample2(d, n):
    assert n == len(d) * 2 or n == len(d) * 2 - 1
    dz = zeros((len(d)-1, 2))
    dz[:,0] = d[0:-1]
    dz[:,1] = 0.5 * d[0:-1] + 0.5 * d[1:]
    #  print dz, ravel(dz)
    if n == len(d) * 2 - 1:
        return array(list(ravel(dz)) + [d[-1]])
    else:
        return array(list(ravel(dz)) + [d[-1],-.5*d[-2] + 1.5*d[-1]])

def upsampleN(d, n):
    l = n
    s = []
    while l != len(d):
        assert l > 0
        s.append(l)
        l = (l+1) / 2
    while s:
        d = upsample2(d, s.pop())
    assert len(d) == n
    return d
    

class scalespace:
    def __init__(self, orig_data, s):
        self.make_scale_space(orig_data, s)
        self.expand_scales()
    def make_scale_space(self, orig_data, s):
        self.scalestep = s
        smoothlimit = 20
    
        d = array(orig_data)
        self.scales = [(0,array(d))]
        octave = -1
        srcimg = 0
        samplelevel = 0
        dsigma = 2**-1
    
        def smoothfactor(s):
            print dsigma
            return sqrt(s**2 - dsigma**2)/2**samplelevel
    
        while 1:
            for i in range(1, s+1):
                targetsigma = 2**(octave + float(i)/float(s))
                while smoothfactor(targetsigma) > smoothlimit:
                    srcimg += 1
                    print "resample %d - %d" % (srcimg,len(self.scales))
    
                    
                    d = downsample(self.scales[srcimg][1])
                    samplelevel = self.scales[srcimg][0] + 1
                    dsigma *= 2**(1.0/s)
                    if len(d) <= 15:
                        return
                    
                sigma = smoothfactor(targetsigma)
                print "Want %.2f, data is already %.2f (%d), so smooth %.2f" % \
                    (targetsigma, dsigma, len(d), sigma)
                self.scales.append((samplelevel, gaussian(d, sigma)))
            
            octave += 1

    def expand_scales(self):
        self.expandscales = []
        for i in range(len(self.scales)):
            self.expandscales.append(self._getidx(i))

    def _getidx(self, si):
        si = int(si)
        if si < 0: si = 0
        elif si >= len(self.scales): si = len(self.scales)-1
        n = len(self.scales[0][1])
        #xassert(all(self.scales[si][1]>0))
        return upsampleN(self.scales[si][1], n)

    def get(self, s):
        getidx = lambda i: self.expandscales[int(i)]
        #getidx = self._getidx
        if s <= 0: s = 0.5 / self.scalestep
        si = (log(s)/log(2) + 1) * self.scalestep

        si1, si2 = floor(si), ceil(si)
        if si1 == si2: return getidx(si1)
        
        d1 = getidx(si1)
        d2 = getidx(si2)
        p = si % 1
        assert 0 < p < 1
        return p * d2 + (1-p) * d1


