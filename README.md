# OPEN CHANNEL HYDRAULICS

## Aim
This is a repository to solve multiple problems related to open channel hydraulics.


## Uniform Flow
Assuming uniform flow in a circular, no circular (rectangular, triangular and trapezoidal) or irregular (e.g. river) channel, the **Manning equation** is used to calculate:

1. Normal water depth

2. Flow discharge

3. Channel slope

4. Manning roughness coefficient

### Structure

#### `chydraulics`

This directory contains:

- `clib.py`: This is a python library with functions to calculate multiple variables related to open channel flow.

- `ufclass.py`: This is a python library with the `UniformClass` to solve problems related to **uniform flow** in regular (e.g. circular, trapezoidal, etc) channels. Here the class `UniformFlow` call functions from `clib.py`. 

- `uficlass.py`: This is a python library with the `UniformClassI` to solve problems related to **uniform flow** in irregular (e.g. rivers) channels. Here the class `UniformFlowI` call functions from `clib.py`. 

#### `test_uf` y `test_ufi`

These directories contain the following:

- `main.py`: This is the python script that call the class `UniformFlow` or `UniformFlowI` in `/chydraulics/ufclass.py` or  `/chydraulics/uficlass.py` which solves the problem.

- `.json` (regular channels): *JSON* is a friendly format to introduce information to scripts. Note that are multiple files with the  extension `.json`, these files contain information for unform flow problems. All the `.json` contain the same structure and information, so that, the scripts are able to indentify which of the four problems need to be solve. The structure of `.json` files is:

   `"ST"`:
    Section type. It can be egual to `1` (Circular cross section) or  `2` (Non circular cross section). [*mandatory*]

   `"US"`:
    Unit measure system. It can be egual to `"IS"` (International system) or  `"BG"` (English system). [*mandatory*]

   `"Q"`:
    Flow discharge. It is given in  `"IS"` or `"BG"`. It is not given for **problem 2**.

   `"So"`:
    Channel slope. It is not given for **problem 3**.

   `"n"`:
    Manning roughness coefficient. It is not given for **problem 4**.

   `"b"`:
    Channel width. It is given in  `"IS"` or `"BG"`. It is  `"b":""` for circular channels.

   `"theta1"`:
    Inclination angle of channel left side with respect to the horizontal. It is given in degrees. It is `"theta1":""` for circular channels.

   `"theta2"`:
    Inclination angle of channel right side with respect to the horizontal. It is given in degrees. It is `"theta2":""` for circular channels.

   `"y"`:
    Normal water depth. It is given in  `"IS"` or `"BG"`. It is not given for **problem 1**.

   `"r"`:
    Radius of a circular channel. It is given in  `"IS"` or `"BG"`. It is `"r":""` for non circular channels.


- `.json` (iregular channels)

   `"US"`:
    Unit measure system. It can be egual to `"IS"` (International system) or  `"BG"` (English system). [*mandatory*]

   `"Q"`:
    Flow discharge. It is given in  `"IS"` or `"BG"`. It is not given for **problem 2**.

   `"So"`:
    Channel slope. It is not given for **problem 3**.

   `"xs"`:
   List of `x` coordinates of the `n` points that shape the cross section.

   `"ys"`:
   List of `y` coordinates of the `n` points that shape the cross section.

   `"ns"`:
   List of Manning roughness coefficients for `n-1` segments (between two point) in the cross section. It is not given for **problem 4**, empty list.

   `"y"`:
   Normal water depth. It is given in  `"IS"` or `"BG"`. It is not given for **problem 1**.

### How to execute it
1. Clone the repo as

  `git clone git@github.com:lamhydro/chydraulics.git`

2. Go into `chydraulics/test_uf` directory:

  `cd chydraulics/test_uf`

3. Execute the code as (e.g.):

  `./main.py uniformFlow_w1020.json`

