# Direct Imaging of Irregular Satellite Disks in Scattered Light

Source code from Nassif and Tamayo (2019): https://arxiv.org/abs/1910.09377. Computes irregular satellite collisional evolution following models from Kennedy and Wyatt (2011). The jupyter\_examples folder provides the notebooks used to generate each of the figures in Nassif and Tamayo (2019), and  demonstrate how to use the API in swarms.py.

## Usage

Two main objects are available in the API, `CollSwarm`, which represents the irregular satellite swarm, as well as `SizeDistribution` which allows for statistical analysis of a swarm's size distribution. `CollSwarm` is the main object that the user should use. A irregular satellite swarm only makes sense in the context of a host planet and its host star, as such we much specify these two objects when creating a new swarm. A star requires its luminosity, mass, temperature, and distance from Earth. Once a star is created, a planet can be specified. A planet requires its host star object, a mass, and a semi-major axis. Additional optional information includes its geometric albedo, its radius, its metalicity, and its age.  
```python
a5 = Star(L = 20 * 3.828e26, M = 2.1 * 2e30, T = 8620, d = 63.4 * 9.4607e15)
jupiter = Planet(star = a5, M = 317.8 * 5.972e24, a = 5.2044 * 1.496e11, Q = 0.5, R = None, 
                    Z = '002', age = 1.e10 * 3.154e7)
``` 
Note that all entries are in SI units. This is true for every function in the API.

Now we can create our swarm,
```python
s = CollSwarm(a5, jupiter, M0 = 0.002 * 7.35e22)
``` 
The collision swarm only requires its initial mass, `M0` to be specified, however there are multiple customization options available for the user. These include, among others, the maximum size of an object in the swarm, the density, stranding corrections, and its age. A full list of its function parameters, as well as the default values, can be found in the docstring of the class in `swarms.py`. Creating a `CollSwarm` object also creates a `SizeDistribution` object attached to the swarm. One can access it by invoking `s.swarm`. 