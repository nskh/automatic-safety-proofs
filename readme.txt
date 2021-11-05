# Abstract for submission
This artifact serves as a complement to our TACAS submission "Automating geometric proofs of collision avoidance via active corners". In our paper, we present a novel approach for efficient, automated collision avoidance verification for convex objects moving in the plan with piecewise trajectories of a single variable. Our approach constructs a formulation of a "safe region", the set of obstacle locations where a collision will not occur given an object and trajectory. An intuitive definition of this safe region is a quantified condition of collision avoidance for the entire trajectory, though the quantifiers in that statement make it inefficient to check at runtime. Our paper proposes a way to construct (provably) equivalent, unquantified representations of such safe regions automatically. The artifact submitted implements our method for the examples handled by our paper. It was used to generate the applications discussed towards the end of our paper, where we apply our method to construct and visualize safe regions corresponding to those in past work by Jeannin et al 2015 and Adler et al 2019. The artifact is a set of Python tools and examples which construct and visualize safe regions given symmetric convex polygons and piecewise planar trajectories that are functions of either x or y. Two types of visualization are supported: 1) approximate, native Python grid visualization of safe and unsafe points (Sympy does not support directly plotting the region corresponding to the boolean safe region formula), or 2) Mathematica plotting of the safe region via Mathematica command output. The examples close by printing Mathematica code that can be copy-pasted to visualize the safe region in Mathematica, if desired.

# README
This README walks through steps required to run the examples for our TACAS submission "Automating geometric proofs of collision avoidance via active corners". Our artifact automatically constructs and visualizes a formulation of a "safe region", the set of obstacle locations where a collision will not occur given an object and trajectory. Commands to run in the terminal begin with a $; do not actually include the "$" character in your terminal command.

On the VM, open a terminal. Download and unzip a zip file of the directory:
TODO(nishant): change this to reflect a different zip file location
$ wget https://github.com/nskh/automatic-safety-proofs/archive/clean_artifact.zip
$ unzip clean_artifact.zip

CD into the directory
$ cd automatic-safety-proofs-clean_artifact

Install dependencies
$ chmod +x install_pip.sh
$ ./install_pip.sh

Four examples are included. The first is a simple toy example to illustrate the functionality of our tool. The second and third are the applications included in our submitted paper. The last example is a piecewise function f(y) that is equivalent to a segment of the second trajectory. When you run the examples, they output three plots (the fourth outputs two plots) and two blocks of text. 

NOTE: in order to continue the example, you must close each plot once it has opened; plots block the continued execution of the example.

The first plot is the polygon object used for unsafe region computation. The second plot is the trajectory of the object center, used to compute the unsafe region. After you close the second plot, a boolean formulation for the unsafe region will be printed to the console. Then there will be a small disclaimer that the next plot may take some time to generate, as in Sympy we cannot plot the boolean formula directly and must repeatedly plot a series of points in a grid to visualize the safe region. The third plot is this grid of points, where blue points are safe and red points are unsafe. Once the third plot is closed, the example concludes by printing Mathematica code that can plot a better visualization of the unsafe region. 

- The first example is a square moving on a piecewise sinusoidal-then-linear path.
- The second example is a rectangle moving on a piecewise parabolic-then-linear path, as in Jeannin 2015.
- The third example is a hexagon moving on a piecewise circular-then-linear path, as in Adler 2019.
- The fourth example is similar to the second, with a parabolic-then-linear path described as a piecewise function f(y). Due to limitations of Sympy, there is no plot of this trajectory in the example, though the computation of the unsafe region works the same way it does in the first three examples.

To run examples:
$ python3 example_1.py
$ python3 example_2.py
$ python3 example_3.py
$ python3 example_4.py