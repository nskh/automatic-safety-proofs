# Automated Geometric Safety Proofs

This README walks through steps required to run the examples for our CAV submission "Automating geometric proofs of collision avoidance via active corners". Our artifact automatically constructs and visualizes a formulation of a "safe region", the set of obstacle locations where a collision will not occur given an object and trajectory. Commands to run in the terminal begin with a $; do not actually include the "$" character in your terminal command.

On the VM, open a terminal. Download and unzip a zip file of the directory:

`$ wget https://github.com/nskh/automatic-safety-proofs/archive/clean_artifact.zip`

`$ unzip clean_artifact.zip`

CD into the directory

`$ cd automatic-safety-proofs-clean_artifact`

Install dependencies

`$ chmod +x install_pip_local.sh`

`$ ./install_pip_local.sh`

Four examples are included. The first is a simple toy example to illustrate the functionality of our tool. The second and third are the applications included in our submitted paper. The last example is a piecewise function f(y) that is equivalent to a segment of the second trajectory. When you run the examples, they output three plots (the fourth outputs two plots) and two blocks of text. 

NOTE: in order to continue the example, you must close each plot once it has opened; plots block the continued execution of the example.

The first plot is the polygon object used for unsafe region computation. The second plot is the trajectory of the object center, used to compute the unsafe region. After you close the second plot, a boolean formulation for the unsafe region will be printed to the console. Then there will be a small disclaimer that the next plot may take some time to generate, as in Sympy we cannot plot the boolean formula directly and must repeatedly plot a series of points in a grid to visualize the safe region. The third plot is this grid of points, where blue points are safe and red points are unsafe. Once the third plot is closed, the example concludes by printing Mathematica code that can plot a better visualization of the unsafe region. 

- The first example is a square moving on a piecewise sinusoidal-then-linear path.
- The second example is a rectangle moving on a piecewise parabolic-then-linear path, as in Jeannin 2015.
- The third example is a hexagon moving on a piecewise circular-then-linear path, as in Adler 2019.
- The fourth example is similar to the second, with a parabolic-then-linear path described as a piecewise function f(y). Due to limitations of Sympy, there is no plot of this trajectory in the example, though the computation of the unsafe region works the same way it does in the first three examples.

To run examples:

`$ python3 example_1.py`

`$ python3 example_2.py`

`$ python3 example_3.py`

`$ python3 example_4.py`