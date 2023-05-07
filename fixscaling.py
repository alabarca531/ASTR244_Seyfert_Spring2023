import matplotlib.pyplot as plt
import math
import numpy as np
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import os
from glob import glob
import warnings
warnings.filterwarnings("ignore")
# warnings.filterwarnings("default")
import matplotlib.gridspec as gridspec

# Real functions that work 
def loadImageData(filename, vbounds = None, graph=True, title = None, save=False):
    if title is None: title = filename
    print(f"Loading {filename}")
    image_data = fits.getdata(filename)
    
    if graph:
        plt.figure()
        if vbounds is not None: plt.imshow(np.log10(image_data), cmap='gray', vmin = vbounds[0], vmax = vbounds[1])
        else: plt.imshow(np.log10(image_data), cmap='gray')
        print('MEAN', np.nanmean(np.log10(image_data)))
        plt.colorbar()
        # plt.grid()
        plt.title(f'{title}')
        if save:
            plt.savefig(f'{title}LoadImage{np.random.randint(500)}.png')

    
    return image_data
    
def distanceFinder(point, center):
    return np.sqrt((center[0]-point[0])**2+(center[1]-point[1])**2)

def plotGalaxy(image, aproxCent = None, vbounds = None, graph = True, radius = 600, title = 'approximate zoom-in', filename=None, save = False):
    image_approximate = image[2200:3200, 1400:2500]
    if aproxCent is not None: 
        xmin, xmax = aproxCent[0]-radius, aproxCent[0]+radius
        ymin, ymax = aproxCent[1]-radius, aproxCent[1]+radius
        image_approximate = image[xmin:xmax, ymin:ymax]
    # center = np.where(image == np.amax(image[-2000:]))
    # print(center)
    # plt.scatter(center[0], center[1], color='yellow')
    # center=(1000,1000)
    # image_smaller = image[int(center[0]):int(center[0])+528, int(center[1])-528:int(center[1])+528]
    if graph:
        plt.figure()
        if vbounds is not None: plt.imshow(np.log10(image_approximate), cmap='gray', vmin = vbounds[0], vmax = vbounds[1])
        else: plt.imshow(np.log10(image_approximate), cmap='gray')
        plt.colorbar()
        plt.title(title)
        plt.scatter(radius, radius, color='red', marker="+")
        if save:
            if filename is not None: plt.savefig(f'{filename}.png')
            

def graph1Dplot(data, axis, title=None, galaxy=None):
    if axis is None: fig, axis = plt.subplots(1,1, figsize = (10,5), sharey=True)
    axis.scatter(data[0], data[1], label = galaxy, alpha = .5)
    axis.legend()
    if title is not None: axis.set_title(title)
    axis.set_xlabel("Distance from brightest point (pixels)")
    axis.set_ylabel("Intensity (Flux)")
    axis.grid(False)

def graph1DSeperate(graphcolor, ax=None, **graphdata):
    if ax is None: fig, ax = plt.subplots(1,1, figsize = (10,5), sharey=True)
    for name, data in graphdata.items():
        ax.scatter(data[0], data[1], label = f'NGC {name[-4:]}', alpha = .5)
    ax.legend()
    ax.set_title(f"Comparison of {graphcolor} Magnitudes as a function of distance")
    ax.set_xlabel("Distance from brightest point (pixels)")
    ax.set_ylabel("Intensity (DN)")
    ax.grid(False)
    

def graphMaker(image, aproxCent, approxHa=None, vbounds = None, graph = True, title = "Approximate Galaxy Zoom-In", radius_gal = 600, radius_core = 50, radius_clust = 40, save = False):
    if approxHa is None: approxTarget, radius_target = aproxCent, radius_core
    else: approxTarget, radius_target = approxHa, radius_clust
     #Crop image down to galaxy size
    xmin, xmax = aproxCent[0]-radius_gal, aproxCent[0]+radius_gal
    ymin, ymax = aproxCent[1]-radius_gal, aproxCent[1]+radius_gal

    image_cropped = image[xmin:xmax, ymin:ymax]
    minval = np.amin(np.log10(image_cropped))
    maxval = np.amax(np.log10(image_cropped))

    if vbounds is not None: minval, maxval = vbounds
    xmin_target, xmax_target = aproxCent[0]-radius_target, aproxCent[0]+radius_target
    ymin_target, ymax_target = aproxCent[1]-radius_target, aproxCent[1]+radius_target
    image_target_cropped = image[xmin_target:xmax_target, ymin_target:ymax_target]
    target_image = image_target_cropped

    minval_target = np.amin(np.log10(target_image))
    maxval_target = np.amax(np.log10(target_image))
    print(minval_target, maxval_target)


    # Now that we have the general idea of where the center is, the max value in the array can be used as the absolute center
    center = np.where(target_image == np.amax(target_image))
    cool_distance = np.array([])
    cool_brightness = np.array([])
    # fig, ax = plt.subplots()
    for i,v in enumerate(target_image):
        for j,u in enumerate(v):
            cool_distance = np.append(cool_distance, distanceFinder(center, (i,j)))
            # print(cool_distance)
            cool_brightness = np.append(cool_brightness, u)
            # cool = [distanceFinder(center, (i,j)), u]
        print(f'{i} finished out of {len(target_image)}', end='\r')
    # print(len(cool_distance), len(cool_brightness))
            # print(f'{cool} given the coordinates were ({i},{j})')
    result = np.array([cool_distance, cool_brightness])
    if graph: 
        fig = plt.figure(tight_layout=True)
        gs = gridspec.GridSpec(2, 2)
        ax_gal = fig.add_subplot(gs[0,0])
        ax_zoom = fig.add_subplot(gs[0,1])
        ax_graph = fig.add_subplot(gs[1,:])
        
        ax_gal.imshow(np.log10(image_cropped), cmap='gray', vmin=minval, vmax=maxval)
        
        ax_gal.scatter(radius_gal, radius_gal, color='red', marker="+")
        # ax_gal.scatter(approxHa[0], approxHa[1], color='green', marker = '.')
        ax_zoom.imshow(np.log10(target_image), cmap='gray', vmin=minval_target, vmax=maxval_target)
        ax_zoom.scatter(center[1], center[0], color='green', marker='o')
        
        graph1Dplot(result, ax_graph)
        graph
        # ax_graph.scatter(cool_distance, cool_brightness, marker='.')
        
        fig.suptitle(title)
    return result


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RangeSlider



# In[14]:


# generate a fake image
def fixScaling2(image1, image2, title=None):
    originalimage1 = image1[0]
    originalimage2 = image2[0]
    scale1 = image1[1]
    scale2 = image2[1]
    img1 = np.log10(originalimage1)
    img2 = np.log10(originalimage2)
    img = np.log10(originalimage1/originalimage2)

    fig, axs = plt.subplots(2, 2, figsize=(10, 5))
    fig.subplots_adjust(bottom=0.25)
    axs[0,0].imshow(img1, vmin = scale1[0], vmax = scale1[1], cmap='gray')
    axs[0,1].imshow(img2, vmin = scale2[0], vmax = scale2[1], cmap='gray')
    im = axs[1,0].imshow(img,cmap='gray' )
    axs[1,1].hist(img.flatten(), bins='auto')
    
    axs[1,1].set_title('Histogram of pixel intensities')
    if title is not None:
        fig.suptitle(f'{title}')

    # Create the RangeSlider
    slider_ax = fig.add_axes([0.20, 0.1, 0.60, 0.03])
    slider = RangeSlider(slider_ax, "Threshold", np.nanmin(img), np.nanmax(img))

    # Create the Vertical lines on the histogram
    lower_limit_line = axs[1,1].axvline(slider.val[0], color='k')
    upper_limit_line = axs[1,1].axvline(slider.val[1], color='k')


    def update(val):
        # The val passed to a callback by the RangeSlider will
        # be a tuple of (min, max)

        # Update the image's colormap
        im.norm.vmin = val[0]
        im.norm.vmax = val[1]

        # Update the position of the vertical lines
        lower_limit_line.set_xdata([val[0], val[0]])
        upper_limit_line.set_xdata([val[1], val[1]])

        # Redraw the figure to ensure it updates
        fig.canvas.draw_idle()


    slider.on_changed(update)
    plt.show()
    
def fixScaling(originalimage, title=None):
    img = np.log10(originalimage)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    fig.subplots_adjust(bottom=0.25)
    im = axs[0].imshow(img,cmap='gray' )
    axs[1].hist(img.flatten(), bins='auto')
    
    axs[1].set_title('Histogram of pixel intensities')
    if title is not None:
        fig.suptitle(f'{title}')

    # Create the RangeSlider
    slider_ax = fig.add_axes([0.20, 0.1, 0.60, 0.03])
    slider = RangeSlider(slider_ax, "Threshold", np.nanmin(img), np.nanmax(img))

    # Create the Vertical lines on the histogram
    lower_limit_line = axs[1].axvline(slider.val[0], color='k')
    upper_limit_line = axs[1].axvline(slider.val[1], color='k')


    def update(val):
        # The val passed to a callback by the RangeSlider will
        # be a tuple of (min, max)

        # Update the image's colormap
        im.norm.vmin = val[0]
        im.norm.vmax = val[1]

        # Update the position of the vertical lines
        lower_limit_line.set_xdata([val[0], val[0]])
        upper_limit_line.set_xdata([val[1], val[1]])

        # Redraw the figure to ensure it updates
        fig.canvas.draw_idle()


    slider.on_changed(update)
    plt.show()
    
# From Notebook

aproxCent_4051 = (1380, 972)
aproxCent_3310 = (1032, 1212)
    
    
    
    
    
path_blue_4051 = 'data/5.7.23_7am_oldimages/Mean NGC4051 LongB.FIT'
path_red_4051 = 'data/5.7.23_7am_oldimages/Mean NGC4051 LongR.FIT'
path_green_4051 = 'data/5.7.23_7am_oldimages/Mean NGC4051 LongG.FIT'
path_ha_4051 = 'data/5.7.23_3am/Mean NGC4051 halpha.FIT'
path_s2_4051 = 'data/5.7.23_3am/Mean NGC4051 sii.FIT'
path_o3_4051 = 'data/5.7.23_3am/Mean NGC4051 oiii.FIT'


path_blue_3310 = 'data/5.7.23_3am/Mean ngc3310.GSC blue.FIT'
path_red_3310 = 'data/5.7.23_3am/Mean ngc3310.GSC red.FIT'
path_green_3310 = 'data/5.7.23_3am/Mean ngc3310 green.GSC.FIT'
path_ha_3310 = 'data/5.7.23_3am/Mean NGC3310 halpha.FIT'

path_s2_3310 = 'data/5.7.23_3am/Mean NGC3310 sii.FIT'
path_o3_3310 = 'data/5.7.23_3am/Mean NGC3310 oiii.FIT'


# data_green_4051 = loadImageData(path_green_4051,vbounds=(10**.4493, 10**.476), graph = False, title='NGC 4051 -- G',save=False)
data_ha_4051 = loadImageData(path_ha_4051,vbounds=(10**.4493, 10**.476), graph = False, title='NGC 4051 -- G',save=False)
# data_green_3310 = loadImageData(path_green_3310,vbounds=(10**.4493, 10**.476), graph = False, title='NGC 3310 -- G',save=False)
data_ha_3310 = loadImageData(path_ha_3310,vbounds=(10**.4493, 10**.476), graph = False, title='NGC 3310 -- G',save=False)


# data_blue_4051 = loadImageData(path_blue_4051, (2,3), graph = False, title='NGC 4051 -- B',save=False)
data_s2_4051 = loadImageData(path_s2_4051, (2,3), graph = False, title='NGC 4051 -- B',save=False)
# data_blue_3310 = loadImageData(path_blue_3310, (2,3), graph = False, title='NGC 3310 -- B',save=False)
data_s2_3310 = loadImageData(path_s2_3310, (2,3), graph = False, title='NGC 3310 -- B',save=False)


# data_red_4051 = loadImageData(path_red_4051,(2,3), graph = False, title='NGC 4051 -- R',save=False)
data_o3_4051 = loadImageData(path_o3_4051,(2,3), graph = False, title='NGC 4051 -- R',save=False)
# data_red_3310 = loadImageData(path_red_3310,(2,3), graph = False, title='NGC 3310 -- R',save=False)
data_o3_3310 = loadImageData(path_o3_3310,(2,3), graph = False, title='NGC 3310 -- R',save=False)


radius_gal_3310 = 150
xmin_3310, xmax_3310 = aproxCent_3310[0]-radius_gal_3310, aproxCent_3310[0]+radius_gal_3310
ymin_3310, ymax_3310 = aproxCent_3310[1]-radius_gal_3310, aproxCent_3310[1]+radius_gal_3310

radius_gal_4051 = 200
xmin_4051, xmax_4051 = aproxCent_4051[0]-radius_gal_4051, aproxCent_4051[0]+radius_gal_4051
ymin_4051, ymax_4051 = aproxCent_4051[1]-radius_gal_4051, aproxCent_4051[1]+radius_gal_4051

# image_cropped = image[xmin:xmax, ymin:ymax]

# fixScaling(data_green_4051[xmin:xmax, ymin:ymax])
# fixScaling(data_ha_4051, title='ha_4051')
ha_4051_scaling = (2.542, 2.620)
# fixScaling(data_green_3310[xmin:xmax, ymin:ymax])
# fixScaling(data_ha_3310, title='ha_3310')
ha_3310_scaling = (1.956, 2.435)

# fixScaling(data_green_4051[xmin:xmax, ymin:ymax])
# fixScaling(data_o3_4051, title='o3_4051')
o3_4051_scaling = (2.546, 2.585)
# fixScaling(data_green_3310[xmin:xmax, ymin:ymax])
# fixScaling(data_o3_3310, title='o3_3310')
o3_3310_scaling = (2.261, 2.485)

# fixScaling(data_green_4051[xmin:xmax, ymin:ymax])
# fixScaling(data_s2_4051, title='s2_4051')
s2_4051_scaling = (2.391, 2.437)
# fixScaling(data_green_3310[xmin:xmax, ymin:ymax])
# fixScaling(data_s2_3310, title='s2_3310')
s2_3310_scaling = (1.929, 2.447)

# fixScaling(data_s2_4051[xmin_4051:xmax_4051, ymin_4051:ymax_4051], 
#            data_ha_4051[xmin_4051:xmax_4051, ymin_4051:ymax_4051], title='s2/ha 4051')
# fixScaling(data_o3_4051[xmin_4051:xmax_4051, ymin_4051:ymax_4051], 
#             data_ha_4051[xmin_4051:xmax_4051, ymin_4051:ymax_4051], title='o3/ha 4051')
# fixScaling(data_s2_3310[xmin_3310:xmax_3310, ymin_3310:ymax_3310], 
#             data_ha_3310[xmin_3310:xmax_3310, ymin_3310:ymax_3310], title='s2/ha 3310')
# fixScaling(data_o3_3310[xmin_3310:xmax_3310, ymin_3310:ymax_3310], 
#             data_ha_3310[xmin_3310:xmax_3310, ymin_3310:ymax_3310], title='o3/ha 3310')

# fixScaling2((data_s2_4051, s2_4051_scaling),
#             (data_ha_4051, ha_4051_scaling),
#             title='s2/ha 4051')
scaling_sh_4051 = (-.214, -.152)
# fixScaling2((data_o3_4051, o3_4051_scaling),
#             (data_ha_4051, ha_4051_scaling),
#             title='o3/ha 4051')
scaling_oh_4051 = (-.100, -0.002)
# fixScaling2((data_s2_3310, s2_3310_scaling),
#             (data_ha_3310, ha_3310_scaling),
#             title='s2/ha 3310')
scaling_sh_3310 = (-0.025, 0.092)
fixScaling2((data_o3_3310[xmin_3310:xmax_3310, ymin_3310:ymax_3310], o3_3310_scaling),
            (data_ha_3310[xmin_3310:xmax_3310, ymin_3310:ymax_3310], ha_3310_scaling),
            title='o3/ha 3310')
scaling_oh_3310 = (-0.0119, 0.2870)
# fixScaling2(data_o3_4051, data_ha_4051, title='o3/ha 4051')
# fixScaling2(data_s2_3310, data_ha_3310, title='s2/ha 3310')
# fixScaling2(data_o3_3310, data_ha_3310, title='o3/ha 3310')