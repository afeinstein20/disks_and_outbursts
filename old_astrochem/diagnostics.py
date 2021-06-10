import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *

COLORMAP = 'viridis'

def one_to_one(ax, low, upp, delta):
    """ Adds a 1 to 1 line to a subplot. """
    ax.plot(np.arange(low, upp, delta),
            np.arange(low, upp, delta),
            c='r', linewidth=3, alpha=0.4)
    return


def time_step_compare(t1, t2, r, zR, locs, 
                      opacity, density,
                      time, UV_lum,
                      UV=1, save=False):
    """
    Creates a diagnostic plot comparing the results
    from two different tables (t1 & t2) at given indices
    (locs). The top plot will be the injected UV flux.
    The middle plots will be the different timesteps.
    The bottom plot will be the initial compared to the 
    final conditions.

    t1 : astropy.table.Table without flare.
    t2 : astropy.table.Table with flare.
    r  : radial distance evaluated at.
    zR : z/r ratio.
    locs : regions to be evaluated.
    time : time array.
    UV_lum : injected UV luminosity array.

    UV : specifies which table to use to plot the flux.
         UV = 0 --> Uses the first input table.
         UV = 1 --> Uses the second input table.
    save : gives the option to save the figure or not.
    """
    cmap = cm.get_cmap(COLORMAP, len(locs)+1)
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3]
        colors.append(matplotlib.colors.rgb2hex(rgb))
    colors = np.array(colors)

    fig = plt.figure(tight_layout=True, figsize=(14,15))
    gs  = gridspec.GridSpec(int(len(locs)/2)+2, 2)
    # Plotting the UV flux
    ax = fig.add_subplot(gs[0,:])

    if UV==0:
        t = t1
    else:
        t = t2

    ax.plot(time, UV_lum, 'k')

    for c, l in enumerate(locs):
        ax.plot(time[l], UV_lum[l], '.',
                c=colors[c], ms=20, label=np.round(time[l],3))

    ax.set_yscale('log')
    ax.set_ylabel('UV Luminosity')
    ax.set_xlabel('Time [Hours]')
    ax.legend(ncol=3)
    ax.set_title('R = {}AU, z/R = {}'.format(r, zR))

    # Sets marker size for plot
    ms = 10

    # Sets the middle plots
    rowlim = np.arange(1, int(len(locs)/2)+1, 1, dtype=int)
    collim = np.arange(0, 2, 1, dtype=int)
    
    l = 0
    for row in rowlim:
        for col in collim:
            ax = fig.add_subplot(gs[row, col])
            one_to_one(ax, 0, 1000, 1)
            
            ax.plot( t1[locs[l-1]], t2[locs[l-1]], '.',
                     c=colors[l-1], ms=ms )
            ax.plot( t1[locs[l]], t2[locs[l]], '.',
                     c=colors[l], ms=ms)

            ax.set_yscale('log')
            ax.set_xscale('log')

            if col == 0:
                ax.set_ylabel('Abundance w/ Flare')

            ax.set_xlim(-19, 1)
            ax.set_ylim(-19, 1)
            
            l += 1


    # Sets the final plot
    ax = fig.add_subplot(gs[int(len(locs)/2)+1,:])
    one_to_one(ax, 0, 1000, 1)

    ax.plot(t1[locs[0]], t1[locs[0]], '.', c=colors[0], ms=ms)
    ax.plot(t1[locs[-1]], t1[locs[-1]], '.', c=colors[-2], ms=ms)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Abundance w/o Flare')
    ax.set_ylabel('Abundance w/ Flare')
    ax.set_xlim(-19, 1)
    ax.set_ylim(-19, 1)

    plt.subplots_adjust()

    if save is False:
        return fig

    else:
        plt.savefig('r{}_zR{}_{}_{}_diagnostic.png'.format(r, zR, opacity, density), 
                    rasterize=True,
                    bbox_inches='tight', dpi=300)
        plt.close()
    
