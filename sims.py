import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon
from scipy.interpolate import interp1d
from descartes import PolygonPatch
from shapely import geometry
import shapely.vectorized

def parse_mask(x):
    ''' Parses a file containing coordinate points for a nucleus outline.
        Input file is a txt file with coordinate points defining the boundaries
        as ordered pairs (i.e., each row has x coord, y coord)
        Output is a 2-d numpy array, where row 1 has x, row 2 had y'''
    with open(x) as f:
        a = f.read()
    b = a.split('\n')
    c = [i.split(', ') for i in b]
    c.pop()
    d = np.zeros((2, len(c)))
    for i in range(len(c)):
        d[0,i] = float(c[i][0])
        d[1,i] = float(c[i][1])
    return d


def parse_roi(x):
    '''Reads in a .txt file with points defining the region of interest
       within the nucleus. Returns a one-D numpy array.'''
    with open(x) as f:
        a = f.read()
    b = a.split(', ')
    d = np.zeros(len(b))
    for i in range(len(b)):
        d[i] = float(b[i])
    return d


def parse_data(data, offset, damage_index):
    '''Reads in a csv with data on accumulation of parp in nucleus over time.
       This assumes that the first several rows of data (time points) are
       before the laser damage occurs.
       x1 is a 2D array with row 1 time points, row 2 parp concentrations, 
       before the laser damage.
       x2 is a 2D array with row 1 time points, row 2 parp concentrations,
       after the laser damage.
       One can use the mean concentration in x1 to do normalize concentrations
       in x2.
       The offset argument is a float, determining at what time the damage
       occurs, and this value is subtracted from all time points in x2 such
       that the time values here (first row) are time since laser damage.
       data is a .csv made by Johannes -- the first five lines are preamble
       stuff which we don't use; the ROI we are tracking is in the 12th column
       (L in an excel spreadsheet) -- this is obviously not generalizable to
       everyone else's data files.
       damage_index (int) defines the length of the two output arrays --
       the entries 1 - damage_index will end up in x1, damage_index - end
       will end up in x2.
       think of it as the row in the csv where the laser damage occurs.'''
    
    with open(data, 'r') as f:
        a = f.read()
    b = a.split('\n')
    # ignore the first 5 entries of b because it the first five lines of
    # the csv are preamble and contain no data
    d = np.zeros((2, len(b)-5))
    for i in range(len(b)-5):
        c = b[i+4].split(',')
        d[0,i] = float(c[0])
        d[1,i] = float(c[11]) # ROI is 12th col of the csv ( hard coded :( )
    # separate off offset for damage time
    d[0,:] -= offset
    # separate pre- and post-damage data into two separate arrays
    x1 = d[:,:damage_index]
    x2 = d[:,damage_index:]
    return x1, x2


def generate_random_points(N, poly):
    '''Randomly (uniformly distributed) placement of N particles within the
       borders of shapely polygon 'poly'. 
       This generates 3N particles and returns only the first N within poly.'''   
    list_of_points = np.zeros((2, N))
    minx,miny,maxx,maxy = poly.bounds
    counter = 0
    x = np.random.uniform(minx, maxx, N*3)
    y = np.random.uniform(miny, maxy, N*3)
    innie = shapely.vectorized.contains(poly, x, y)
    return x[innie][:N], y[innie][:N]

def D_2_x(D, h):
    ''' Wrapper for transforming the diffusion coefficient into delta x.
        Here x = dx, h = dt (just like in calculus class)
        D = (1/2) * (dx)**2 / (dt) <-> dx = sqrt(2 * D * dt)'''
    x = np.sqrt((2*D*h))
    return x

def update_positions(x_cord, y_cord, mu, sigma, nucleus, roi):
    ''' Function updating the positions of particles along a lattice.
        This assumes the following:
        - step size is distributed according to a Gaussian (mu, sigma)
        - particles can move in x and y (or both), i.e., can move in one
        of eight directions (left, left and up, up, right and up, ...) but
        must move in at least one direction (direction of movement is 50/50)
        - once particles enter the ROI, they are STUCK in the ROI and
        no longer move
        We took the following liberty while modeling:
        - When particles become stuck in the ROI, we no longer store their
        positions; we instead return the number of particles which became
        stuck in the ROI in that time step (this means the number of
        coordinate pairs returned may be smaller than the number passed in)
        - When we simulate movement that moves particles outside the ROI,
        then we do not update the particle's position. This can be thought of
        as similar to reflecting boundary conditions.
        INPUT PARAMS:
        x_cord, y_cord are x, y coordinates to be updated
        mu, sigma define a Gaussian for step size.
        nucleus is a shapely polygon defining the nucleus
        roi is a shapely polygon defining the ROI
        OUTPUT PARAMS:
        x_new[out_roi], y_new[out_roi] are update x, y coords of particles
        which are not in the ROI
        N_stuck is the number (integer) of particles which got stuck in
        the ROI in this time step
        '''
    # Number of particles to simulate moving
    l = len(x_cord)
    # Draws from Gaussian to define step length in x and y (independent)
    x = np.random.normal(mu, sigma, l)
    y = np.random.normal(mu, sigma, l)
    # Draws from uniform distribution to define direction of movement
    # (this has 50% chance of generating -1, 50% of generatin 1)
    # note x and y directions of movement are independent of each other
    fx = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    fy = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    # combine step length and direction
    x *= fx
    y *= fy
    # create updated positions for each particle
    x_new = x_cord + x
    y_new = y_cord + y
    # check to see which of these new updated positions is outside the nucleus
    # (out_nuc is a boolean mask)
    out_nuc = shapely.vectorized.contains(nucleus, x_new, y_new) == False
    # for positions which are not in the nucleus, change updated positions
    # to be the pre-motion positions (i.e., don't change position if the
    # new position would be outside the nucleus)
    x_new[out_nuc] = x_cord[out_nuc]
    y_new[out_nuc] = y_cord[out_nuc]
    # Determine number of particles which have entered the ROI
    # (note -- no coordinates of particles inside the ROI should have been
    # passed into this function, so all particles in the ROI are new to the
    # ROI in this time step)
    in_roi = shapely.vectorized.contains(roi, x_new, y_new)
    N_stuck = np.sum(in_roi)
    # Give boolean mask for particles outside the ROI
    # (so that we can return only positions of these particles)
    out_roi = in_roi == False

    return x_new[out_roi], y_new[out_roi], N_stuck


def init_sim(N, nuc):
    '''A wrapper function for initializing a simulation (x and y init. coords.)
       We found it useful to initialize once then run many simulations,
       as this keeps the amount of initial particle in ROI constant,
       which removes artificial variance when normalizing results.
       We then recommend optimizing and estimating parameters over several
       initial conditions.
       INPUT PARAMS:
       N - number of particles
       nuc - shapely polygon of nucleus
       OUTPUT PARAMS:
       x, y - arrays of positions of N particles within nucleus'''
    x, y = generate_random_points(N, nuc)
    return x, y

def simulate(D, f_mobile, f_bleached, nuc, roi, runtime, x0, y0):
    ''' Function to simulate particle movement for specified runtime.
        This takes in already-initialized points and, then updates positions
        the specified number of itme steps.
        INPUT PARAMS:
        D - diffusion coefficient in microns sq. per second (input D to
        make a D, F parameter sweep easier)
        f_mobile - mobile fraction, the proportion of particles in cell which
        can move
        f_bleached - proportion of particles in ROI which become photobleached
        nuc - shapely polygon object defining nucleus boundaries
        roi - shapely polygon object defining ROI boundaries
        runtime - number of timesteps to update particle positions
        (** I think this is currently in timesteps -- shouldn't it be in secs?)
        x0, y0 - initial positions (should be passed in from init_sim())
        OUTPUT PARAMS:
        all_stuck - (array) time series with counts of particles in ROI
        '''
        
    # length of each time step (in seconds)
    # we hard-coded this in to roughly match Johannes's experiments.
    h = 0.18
    # conversion factor for microns to pixels
    microns_2_pixels = .08677
    # number of particles to simulate 
    # (hard-coded for now, maybe should later be added as argument for
    # larger nucleus sizes)
    N = 12000
    # calculating dx (in pixels) from D (in microns)
    dx = D_2_x(D, h) / microns_2_pixels
    # in_roi is a boolean mask for each particle, whether or not in ROI
    in_roi = shapely.vectorized.contains(roi, x0, y0)
    # out_roi is the compliment for in_roi (also a boolean mask)
    out_roi = in_roi == False
    # N0_roi is the total number of particles in the ROI at time 0
    N0_roi = np.sum(in_roi)
    # stuck is the number of particles in the ROI which did not get bleached
    stuck = int(N0_roi * (1 - f_bleached))
    # initialize all_stuck, a numpy array which will hold the number of
    # particles stuck in the ROI in each time step
    all_stuck = np.zeros(runtime+1)
    all_stuck[0] = stuck
    # N_sim is the number of particles to simulate moving,
    # x, y are the coordinates within the nucleus of these particles
    # (i.e., the proportion of the particles not in the ROI which are mobile)
    N_sim = int((N - N0_roi) * f_mobile)
    x = x0[out_roi][:N_sim]
    y = y0[out_roi][:N_sim]
    
    # iterate particle movements over time
    # in each time step, store the number of particles moved
    for i in range(1,runtime+1):
        # update positions of each particle (n.b. sd of the Gaussian is fixed)
        x, y, N_stuck = update_positions(x, y, dx, 0.001, nuc, roi)
        # N_stuck is the number of particles which got stuck in the ROI
        # in this time step
        # stuck is the overall number of particles stuck in the ROI
        stuck += N_stuck
        # store number of particles stuck in ROI in this time step
        all_stuck[i] = stuck

    return all_stuck


def compute_error(data, data_norm, stuck_time, stuck_norm):
    times = np.zeros(len(data[0,:]))
    y_ = np.zeros(len(times))
    error = np.zeros(len(times))
    for i in range(len(times)):
        dx = (np.abs(stuck_time - data[0,i])).argmin()
        times[i] = stuck_time[dx]
        y_[i] = stuck_norm[dx]
        error[i] = (stuck_norm[dx] - data_norm[i])**2.
    return error #np.sqrt(error)

def compute_Rsq(data, data_norm, stuck_time, stuck_norm, plot_name = None):
    '''Calculates R-squared 
    (technically not correct Rsq due to non-linearity of our curve but w/e)
    INPUT PARAMS:
    data - raw observed data, 2d numpy array, where first row contains time (s)
    data_norm - normalized observed concentration (within ROI) data
    stuck_time - array with time (s) of simulation output
    stuck_norm - array with concentrations (normalized) in ROI from simulation
    plot_name - name of plot (if None, do not plot results)
    OUTPUT PARAMS:
    R squared = 1 - (sum of sq. resid. error) / (sum of sq. total error)'''
    
    # interpf is an interpolation object (from scipy)
    # used to interpolate simulated concentration in ROI for the time points
    # in observed data
    interpf = interp1d(stuck_time, stuck_norm)
    # "predicted" concentrations by simulation at each time step in obs. data
    predict = interpf(data[0,:])
    # residuals
    residus = predict - data_norm
    
    # summed square error (total variance in observed data)
    sse = np.sum((data_norm - np.mean(data_norm))**2)
    # summed square residuals (error in simulation)
    ssr = np.sum((residus)**2)
    
    if not plot_name:
        # plot residuals
        plt.plot(data[0,:], residus, '.')
        # plot zero line for comparison
        plt.plot(np.linspace(0, max(data[0,:]), len(data[0,:])), np.zeros(len(data[0,:])),'-')
        plt.savefig(plot_name)
    
    return 1 - (ssr/sse)