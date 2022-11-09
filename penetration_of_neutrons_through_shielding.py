#!/usr/bin/env python
# coding: utf-8

# # Project 3: Monte Carlo Techniques: Penetration of Neutrons Through Shielding
#     April-May 2022
#     University of Manchester
#     Arnau Duran Mayol
# 
# 
# ### Introduction
# The task in this project is to develop a simulation of penetration of neutrons through a slab of shielding of thickness T, considering only thermal neutrons, and the processes of absorption and scattering.
# The simulation will be run using random number generators and isotropic step generators for three materials (water, lead and graphite).

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as scp
from matplotlib import cm


# ## Random number
# Different methods of generating random numbers are tested before starting the project. Initially, $np.random.uniform()$ is used and it is visualised in a 3D plot to look for the spectral problem.

# In[2]:


def random_num(lim1, lim2, n):
    '''
    Function that generates random numbers
    Args:
        lim1: start value for random number generation
        lim2: end limit for random number generation
        n: number of samples
    Returns:
        random 3D array
    '''
    x_values = np.random.uniform(lim1, lim2, size=n)
    y_values = np.random.uniform(lim1, lim2, size=n)
    z_values = np.random.uniform(lim1, lim2, size=n)
    return x_values, y_values, z_values

x_values, y_values, z_values = random_num(1, 10, 10000)
plt.rcParams["figure.figsize"] = (6,4)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x_values, y_values, z_values, c='g')
plt.show()


# ## Spectral problem
# The majority of random number generator are made of the linear congruential generatirs(LCGs). When LCGs are plotted in 2 or more dimensions, lines or planes can be seen, on which all possible outputs can be found. This is called spectral problem.
# 
# As we can see, this random or pseudorandom generator doesn't present the spectral problem as there are no paralel areas where the points are generated.
#  
# ## RANDSSP generator
# Here another random number generator method is verified, the RANDSSP. This generator presents the spectral problem as we can see clearly the paralel planes.
# 

# In[3]:


get_ipython().run_line_magic('matplotlib', 'notebook')

def randssp(p,q):
    # RANDSSP Multiplicative congruential uniform random number generator.
    # Generates random matrix
    global m, c, x
        
    try: a
    except NameError:
        m = pow(2, 31)
        a = pow(2, 16) + 3
        c = 0
        x = 123456789
    
    try: p
    except NameError:
        p = 1
    try: q
    except NameError:
        q = p
    
    r = np.zeros([p,q])

    for l in range (0, q):
        for k in range (0, p):
            x = np.mod(a*x + c, m)
            r[k, l] = x/m
    
    return r
random_array = randssp(3,10000)
plt.rcParams["figure.figsize"] = (6,4)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(random_array[0, :], random_array[1, :], random_array[2, :], c='b')
plt.show()


# ## Penetration process

# A function is coded to find the number of absorbed neutrons using
# $$n = \frac{\rho N_A}{M}$$
# Where $\rho$ is the density, $M$ the molar mass of the material and $N_A$ the Avogadro's Number.
# Then the mean free path of the neutron is defined as
# $$\lambda=\frac{1}{n\sigma _a}$$
# Where $\sigma _a$ is the cross section and $n$ the number of absorbed neutrons.
# The final output function is
# $$s_i = -\lambda log（a_i)$$
# Where $a_i$ is a random number generated form 0 to 1 and $\lambda$ the mean free path. We write a random number generator that generates samples distributed according to an exponential function.
# 

# In[4]:


d = 1.00 # water density g/cm3
M = 18.01528 # Molar mass g/mol
def absorbing(d, M):
    '''
    Function to find the number of absorbed neutron
    Args:
        d: float for density
        M: float for molar mass
    Returns:
        float for number of absorbed neutrons
    '''
    absorbed_neutrons = d*scp.N_A / M
    return absorbed_neutrons

sigma_a = 0.6652*10**-24 #change from barns to cm2
sigma_s = 103*1.0*10**-24 #change from barns to cm2, there is abscence for scattering

def mfp(d, M, sigma):
    '''
    Function to find the mean free path
    Args:
        d: float for density
        M: float for molar mass
        sigma: float for cross section
    Returns:
        float for mean free path
    '''
    mfp = 1/(absorbing(d, M)*sigma)
    return mfp

def exponential_gen(n, d, M, sigma):
    '''
    Function to find the generator that generates samples distributed according to an exponential function
    Args:
        n: integer for number of samples
        d: float for density
        M: float for molar mass
        sigma: float for cross section
    Returns:
        array of random numbers according to an exponential function
    '''
    a = np.random.uniform(0, 1, n)
    return -mfp(d, M, sigma) * np.log(a)

plt.rcParams["figure.figsize"] = (6,4)
fig = plt.figure()
plt.hist(exponential_gen(10000, d, M, sigma_a), 40) #plotting histogram of 40 bins
plt.xlabel('thickness of water (cm)')
plt.ylabel('number of samples')
plt.show()


# In[5]:


n=1000
Nr, b = np.histogram(exponential_gen(n, d, M, sigma_a), 40)
b_final = (b[:-1]+b[1:])/2
Nr = Nr[np.where(Nr > 0)]
b_final = b_final[np.where(Nr > 0)]
b_error = (b[1:]-b[0:len(b)-1])/2
b_error = b_error[np.where(Nr > 0)]
(coef, covr) = np.polyfit(b_final, np.log(Nr), 1, cov=True)
error_fit = np.sqrt(covr[0][0])

fit = np.polyval(coef, b_final)
#error1 = np.sqrt(covr[0][0]
plt.rcParams["figure.figsize"] = (6,4)
fig = plt.figure()
plt.errorbar(b_final, np.log(Nr), xerr=b_error, fmt='.')
plt.plot(b_final, fit, '--', label='Linear fit')
plt.xlabel('thickness of water (cm)')
plt.ylabel('log(number of samples)')
plt.legend()
plt.show()
print('The gradient value is {:04.4f} +/- {:04.4f} and the  mean free path is {:04.3f} cm.'.format(coef[0], error_fit, mfp(d, M, sigma_a)))


# ## 3D Randomly generated sphere
# What follows next is to generate random arrays of isotropic 3-D unit vectors, using a spherical coordinate system.

# In[6]:


n = 1000

def sphere_vector(n):
    '''
    Function that creates a random sphere
    Args:
        n: size of the sample
    Returns:
        3D array of random values
    '''
    u = np.random.uniform(0, 1, n)
    theta = np.arccos(1- 2*(u))
    phi = np.random.uniform(0, 2*np.pi, n)
    r = mfp(d, M, sigma_a)
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z

#scatter plot to show that the points of the sphere are determined randomly
get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure()
ax = Axes3D(fig)
x,y,z=sphere_vector(n)
ax.scatter(x, y, z, c='black')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()


# As we can see from the  plot, we have obtained a succesfull 3D randomly uniform sphere.

# ## Step function
# 
# Finally, isotropic steps can be generated given that the length $l$ follows an exponential distribution.
# 

# In[7]:


l = np.exp(-1/mfp(d, M, sigma_a + sigma_s))
n=1000

def step_function(l, n):
    '''
    Function that generates 3D random steps
    Args:
        l: attenuation length
        n: number of samples
    Returns:
        3D array of random values
    '''
    u = np.random.uniform(0, 1, n)
    theta = theta = np.arccos(1-2*np.random.uniform(0, 1, n))
    phi = np.random.uniform(0, 2*np.pi, n)
    r = -l * np.log(np.random.uniform(size=n))
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

array = np.empty((n, 3))
j=0
while j<n:
    array[j] = step_function(l, 1)
    j+=1 #adds steps together
x_step = array[:, 0]
y_step = array[:, 1]
z_step = array[:, 2]
get_ipython().run_line_magic('matplotlib', 'notebook')
plt.rcParams["figure.figsize"] = (6,4)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x_step, y_step, z_step, c='black')
ax.set_title('3D plot of 1000 random steps of neutrons in water')
ax.set_xlabel('x(cm)')
ax.set_ylabel('y(cm)')
ax.set_zlabel('z(cm)')
plt.show()


# ## COMPLETE SIMULATION
# A complete simulation of a neutron in different mediums (lead, graphite and water) is coded. Initially a thickness of 10 cm is used for all mediums.
# When we have multiple processes inside the material, the mean free path now becomes
# $$\lambda=\frac{1}{n(\sigma _a + \sigma _s)}$$
# The absorbtion probability to filter and determine the state if each neutron is:
# $$ P_a = \frac {\sigma_a} {\sigma_a + \sigma_s}$$
# 

# In[104]:


n = 1000
thickness = 10 #in cm, initially we use this thickness for all materials

def stop_function(x_step, thickness, u, sigma_a, sigma_s):
    '''
    Function that stops the porcess when a neutron is transmitted, scattered or absorbed
    Args:
        x_step: random array in the x step
        thickness: integer
        u: random number
        sigma_a: absorbtion cross section
        sigma_s: scattering cross section
    Returns:
        Transmission, scattering or absorbtion
    '''
    if x_step > thickness:
        return 1 #transmission
    if x_step < 0:
        return 2 #scattering
    if u < sigma_a / (sigma_a + sigma_s):
        return 3 #absorbtion
    else:
        return 0

def move_function(thickness, sigma_a, sigma_s):
    '''
    Function that moves randomly in a 3D medium
    Args:
        thickness: integer
        sigma_a: absorbtion cross section
        sigma_s: scattering cross section
    Returns:
        3D array of random values
        result: final stage of the neutron
    '''
    sigma = sigma_a + sigma_s
    x_move = [0]
    y_move = [0]
    z_move = [0]
    #first step
    x_step = exponential_gen(1, d, M, sigma)
    y_step = 0
    z_step = 0
    num = np.random.uniform()
    # while loop which adds to the initial step the step function
    while stop_function(x_step, thickness, num, sigma_a, sigma_s) == 0:
        x_move = np.append(x_move, x_step)
        y_move = np.append(y_move, y_step)
        z_move = np.append(z_move, z_step)
        x_vector, y_vector, z_vector = step_function(l, 1)
        x_step = x_step + x_vector
        y_step = y_step + y_vector
        z_step = z_step + z_vector
    x_move = np.append(x_move, x_step)
    y_move = np.append(y_move, y_step)
    z_move = np.append(z_move, z_step)  
    #chosing of final stage of neutron
    if stop_function(x_step, thickness, num, sigma_a, sigma_s) == 1:
        result = 'This process goes through transmission'
    if stop_function(x_step, thickness, num, sigma_a, sigma_s) == 2:
        result = 'This process goes through scattering'
    if stop_function(x_step, thickness, num, sigma_a, sigma_s) == 3:
        result = 'This process goes through absorbtion'
    return (x_move, y_move, z_move, result)

def neutron_class(sigma, n, thickness):
    '''
    Function that finds the quantity of each possible process (transmission, scattering and absorbtion)
    Args:
        thickness: integer
        n: number of samples
        sigma: cross section
    Returns:
        quantity of each possible process
    '''
    transmission = 0
    scattering = 0
    absorbtion = 0
    a = np.random.uniform(size=n)
    new_x = exponential_gen(n, d, M, sigma)
    while len(new_x) > 0:
        u = np.random.uniform(size=len(new_x))
        transmission = transmission + np.count_nonzero(new_x>thickness)
        scattering = scattering + np.count_nonzero(new_x<0)
        absorbtion = absorbtion + np.count_nonzero(u[np.argwhere((new_x>0)&(new_x<thickness))]<(sigma_a/(sigma_a+sigma_s)))
        new_x = np.delete(new_x, np.argwhere((new_x<0) | (new_x>thickness) | (u<sigma_a/(sigma_a+sigma_s))))
        new_x = new_x + step_function(l, len(new_x))[0]
    return transmission, scattering, absorbtion

def percentage_neutrons(transmission, scattering, absorbtion):
    '''
    Function that finds the percentage of transmission, scattering and absorbtion
    Args:
        tranmsission: integer
        scattering: integer
        absorbtion: integer
    Returns:
        percentage of each process
    '''
    total = transmission + scattering + absorbtion
    percentage_t = transmission/total*100
    percentage_s = scattering/total*100
    percentage_a = absorbtion/total*100
    return percentage_t, percentage_s, percentage_a


# ## WATER

# In[97]:


d = 1 # water density g/cm3
M = 18.01528 # Molar mass g/mol
sigma_a = 0.6652*1.0*10**-24 #change from barns to cm2
sigma_s = 103*1.0*10**-24 #change from barns to cm2, there is abscence for scattering
l = np.exp(-1/mfp(d, M, sigma_a+sigma_s))
x_move, y_move, z_move, result = move_function(thickness, sigma_a, sigma_s)
transmission, scattering, absorbtion = neutron_class(sigma_a+sigma_s, n, thickness)
percentage_t, percentage_s, percentage_a = percentage_neutrons(transmission, scattering, absorbtion)

fig = plt.figure(figsize=(9, 4))
ax1 = fig.add_subplot(121, projection='3d')
ax1.scatter(x_move[0], y_move[0], z_move[0], c='blue', label='Starting value')
ax1.scatter(x_move[-1], y_move[-1], z_move[-1], c='red', label='Final value')
ax1.plot(x_move, y_move, z_move)
ax1.set_title('Path of a neutron in water medium')
ax1.set_xlabel('x(cm)')
ax1.set_ylabel('y(cm)')
ax1.set_zlabel('z(cm)')
ax1.legend()

print(result)
print('The mean free path in a system of water is {0:.3f} and the attenuation length {1:.3f}'.format(mfp(d, M, sigma_a+sigma_s), l))
print('The percentage of transmission is {0:.3f}%, of scattering {1:.3f}% and of absorbtion {2:.3f}%'.format(percentage_t, percentage_s, percentage_a))

neutrons = [transmission, scattering, absorbtion]
state = ['Transmission', 'Scattering', 'Absorbtion']
ax2 = fig.add_subplot(122)
ax2.set_title('Distribution of neutrons in water')
ax2.pie(neutrons, labels = state)
plt.show()


# ## LEAD

# In[105]:


d = 11.35  # lead density g/cm3
M = 207.2 # Molar mass g/mol
sigma_a = 0.158*1.0*10**-24 #change from barns to cm2
sigma_s = 11.221 *1.0*10**-24 #change from barns to cm2, there is abscence for scattering
sigma = sigma_a + sigma_s
l = np.exp(-1/mfp(d, M, sigma))
x_move, y_move, z_move, result = move_function(thickness, sigma_a, sigma_s)
transmission, scattering, absorbtion = neutron_class(sigma_a+sigma_s, n, thickness)
percentage_t, percentage_s, percentage_a = percentage_neutrons(transmission, scattering, absorbtion)

fig = plt.figure(figsize=(9, 4))
ax1 = fig.add_subplot(121, projection='3d')
ax1.scatter(x_move[0], y_move[0], z_move[0], c='blue', label='Starting value')
ax1.scatter(x_move[-1], y_move[-1], z_move[-1], c='red', label='Final value')
ax1.plot(x_move, y_move, z_move)
ax1.set_title('Path of a neutron in lead medium')
ax1.set_xlabel('x(cm)')
ax1.set_ylabel('y(cm)')
ax1.set_zlabel('z(cm)')
ax1.legend()

print(result)
print('The mean free path in a system of lead is {0:.3f} cm and the attenuation length {1:.3f} cm'.format(mfp(d, M, sigma_a+sigma_s), l))
print('The percentage of transmission is {0:.3f}%, of scattering {1:.3f}% and of absorbtion {2:.3f}%'.format(percentage_t, percentage_s, percentage_a))

neutrons = [transmission, scattering, absorbtion]
state = ['Transmission', 'Scattering', 'Absorbtion']
ax2 = fig.add_subplot(122)
ax2.set_title('Distribution of neutrons in lead')
ax2.pie(neutrons, labels = state)
plt.show()


# ## GRAPHITE

# In[99]:


d = 1.67  # water density g/cm3
M =  12.01# Molar mass g/mol
sigma_a = 0.0045*1.0*10**-24 #change from barns to cm2
sigma_s = 4.74 *1.0*10**-24 #change from barns to cm2, there is abscence for scattering
l = np.exp(-1/mfp(d, M, sigma_a+sigma_s))
x_move, y_move, z_move, result = move_function(thickness, sigma_a, sigma_s)
transmission, scattering, absorbtion = neutron_class(sigma_a+sigma_s, n, thickness)
percentage_t, percentage_s, percentage_a = percentage_neutrons(transmission, scattering, absorbtion)

fig = plt.figure(figsize=(9, 4))
ax1 = fig.add_subplot(121, projection='3d')
ax1.scatter(x_move[0], y_move[0], z_move[0], c='blue', label='Starting value')
ax1.scatter(x_move[-1], y_move[-1], z_move[-1], c='red', label='Final value')
ax1.plot(x_move, y_move, z_move, c='green')
ax1.set_title('Path of a neutron in graphite medium')
ax1.set_xlabel('x(cm)')
ax1.set_ylabel('y(cm)')
ax1.set_zlabel('z(cm)')
ax1.legend()

print(result)
print('The mean free path in a system of graphite is {0:.3f} and the attenuation length {1:.3f}'.format(mfp(d, M, sigma_a+sigma_s), l))
print('The percentage of transmission is {0:.3f}%, of scattering {1:.3f}% and of absorbtion {2:.3f}%'.format(percentage_t, percentage_s, percentage_a))

neutrons = [transmission, scattering, absorbtion]
state = ['Transmission', 'Scattering', 'Absorbtion']
ax2 = fig.add_subplot(122)
ax2.set_title('Distribution of neutrons in graphite')
ax2.pie(neutrons, labels = state)
plt.show()


# As we can see from this plots, we have obtained a successfull simulation for neutrons in different mediums, with different distributions of transmission, abosrbtion and scattering.

# ## Variating slab thickness
# To continue with the project, this time we are going to see what the change in thickness does to transmission, absobtion and scattering rates.

# In[100]:


thickness = np.geomspace(0.001,10, 30) #array of thickness values between 0.001 and 10
r_transmission = np.empty(len(thickness))
r_scattering = np.empty(len(thickness))
r_absorbtion = np.empty(len(thickness))

def thickness_variation(sigma, n, thickness):
    '''
    Function that finds transmission, scattering and absorbtion values when thickness changes
    Args:
        sigma: float of cross section
        n: integer for number of samples
        thickness: array of floats
    Returns:
         r_transmission, r_scattering, r_absorbtion: array of floats
    '''
    for i in range(len(thickness)):
        r_transmission[i],_,_ =  neutron_class(sigma, n, thickness[i])
    for i in range(len(thickness)):
        _,r_scattering[i],_ =  neutron_class(sigma, n, thickness[i])
    for i in range(len(thickness)):
        _,_,r_absorbtion[i] =  neutron_class(sigma, n, thickness[i])
    return r_transmission, r_scattering, r_absorbtion

def plot_variation(r_transmission, r_scattering, r_absorbtion, title):
    '''
    Plotting function
    Args:
        r_transmission: array of transmission values
        r_scattering: array of scattering values
        r_absorbtion: array of absorbtion values
        title: title of the plot
    Returns:
         plots
    '''
    plt.rcParams["figure.figsize"] = (14,4)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.suptitle(title)
    ax1.scatter(thickness, r_transmission)
    ax1.set_xlabel('Thickness(cm)')
    ax1.set_ylabel('Rate of transmissions')
    ax2.scatter(thickness, r_scattering)
    ax2.set_xlabel('Thickness(cm)')
    ax2.set_ylabel('Rate of scattering')
    ax3.scatter(thickness, r_absorbtion)
    ax3.set_xlabel('Thickness(cm)')
    ax3.set_ylabel('Rate of absorbtion')
    plt.show()
    return 0

sigma_a = 0.6652*1.0*10**-24 #change from barns to cm2
sigma_s = 103*1.0*10**-24 #change from barns to cm2, there is abscence for scattering
r_transmission, r_scattering, r_absorbtion = thickness_variation(sigma_a+sigma_s, n, thickness)
water_variation = plot_variation(r_transmission, r_scattering, r_absorbtion, 'Variations of thickness for water')

sigma_a = 0.158*1.0*10**-24 #change from barns to cm2
sigma_s = 11.221 *1.0*10**-24 #change from barns to cm2, there is abscence for scattering
r_transmission, r_scattering, r_absorbtion = thickness_variation(sigma_a+sigma_s, n, thickness)
lead_variation = plot_variation(r_transmission, r_scattering, r_absorbtion, 'Variations of thickness for lead')

sigma_a = 0.0045*1.0*10**-24 #change from barns to cm2
sigma_s = 4.74 *1.0*10**-24 #change from barns to cm2, there is abscence for scattering
r_transmission, r_scattering, r_absorbtion = thickness_variation(sigma_a+sigma_s, n, thickness)
graphite_variation = plot_variation(r_transmission, r_scattering, r_absorbtion, 'Variations of thickness for graphite')


# As we can see from the plots, rate of transmission decreases as thickness increases. On the other hand scattering and absorbtion rates seem to increase with thickness.
# 
# A best fit for transmission rates against thickness in order to obtain the attenuation length is made. Where attenuation length is
# $$l=-\frac {1}{m}$$
# Where $l$ is the attenuation length and $m$ is the gradient.

# In[101]:


def best_fit(r_transmission, title):
    '''
    Function that finds the best fit and plots it
    Args:
        r_transmission: array of transmission values
        title: title of the plot
    Returns:
         plot
         print of attenuation length values
    '''
    plt.rcParams["figure.figsize"] = (5,4)
    fig = plt.figure()
    (coef, covr) = np.polyfit(thickness, np.log(r_transmission), 1, cov=True)
    fit = np.polyval(coef, thickness)
    plt.scatter(thickness, np.log(r_transmission))
    plt.plot(thickness, fit, '--', color='red')
    plt.title(title)
    plt.xlabel('Thickness(cm)')
    plt.ylabel('Rate of transmissions')
    plt.show()
    error_grad = np.sqrt(covr[0][0])
    gradient = coef[0]
    attenuation = -1/gradient

    print('The attenuation length is {0:.3f}  +/- {1:.3f} cm.'.format(attenuation, error_grad))

    return 0

sigma_a = 0.6652*1.0*10**-24 #change from barns to cm2
sigma_s = 103*1.0*10**-24 #change from barns to cm2, there is abscence for scattering
r_transmission, r_scattering, r_absorbtion = thickness_variation(sigma_a+sigma_s, n, thickness)
water_fit = best_fit(r_transmission, 'Tranmission fit for water')
sigma_a = 0.158*1.0*10**-24 #change from barns to cm2
sigma_s = 11.221 *1.0*10**-24 #change from barns to cm2, there is abscence for scattering
r_transmission, r_scattering, r_absorbtion = thickness_variation(sigma_a+sigma_s, n, thickness)
lead_fit = best_fit(r_transmission, 'Tranmission fit for lead')
sigma_a = 0.0045*1.0*10**-24 #change from barns to cm2
sigma_s = 4.74 *1.0*10**-24 #change from barns to cm2, there is abscence for scattering
r_transmission, r_scattering, r_absorbtion = thickness_variation(sigma_a+sigma_s, n, thickness)
graphite_fit = best_fit(r_transmission, 'Tranmission fit for graphite')


# ## BINOMIAL ERRORS
# The code was unsuccesful with error analysis.

# In[80]:


thickness = np.geomspace(0.001,10, 30) #array of thickness values between 0.001 and 10

data = np.empty([10, 3])
def error_samples(sigma_a, sigma_s, n, thickness):
    error_sample = np.empty(())
    for i in range(10):
        data[k, :] =  thickness_variation(sigma, n, thickness)
        std_samples[i],_,_ = np.std(data)
        error_samples = std_samples / mean_samples
    return error_samples
print(error_samples(sigma_a, sigma_s, n, thickness))


# ## Conclusion 
# In conclusion, this project has allowed us to work with random number generators, which are the basis for Monte Carlo method. We have also succesfully coded a random simulation for neutrons going into a medium. These are transmitted, scattered and absorbed in different rates, depending on the thickness and the medium.
# The physics knowledge adquired with this project is a better understanding of scattering of neutrons, mean free path and attenuation lengths for different materials.

# In[ ]:




