# Littlewood Polynomials

A polynomial 

<img src="https://render.githubusercontent.com/render/math?math=\hspace{11cm}\displaystyle{p(x) = \sum_{i = 0}^n a_i x^i,}">

is a *Littlewood polynomial* if all the <img src="https://render.githubusercontent.com/render/math?math=a_i = \pm 1">. The objective of this project is to 
carry out an exploratory analysis with continuous perturbations of Littlewood polynomials, i.e., to see what happens in the neighborhood of these polynomials.
A (Jupyter) notebook is available for anyone interested in conducting their own experiments. This notebook generates a sequence of images which should then be 
converted to video. Each frame is a heatmap of zeros of a certain set of polynomials, where darker color means less zeros in the region, and brighter colors 
means more zeros.

We illustrate what this projects does with a simple example. In [this video](https://www.youtube.com/watch?v=wZZqCccU0wk) we see what happens when the independent 
coefficient goes to infinity. More precisely, for each <img src="https://render.githubusercontent.com/render/math?math=t \in [0, +\infty)"> we have a frame of the 
video, which shows the set of zeros of all perturbed Littlewood polynomials 

<img src="https://render.githubusercontent.com/render/math?math=\hspace{10cm}p(x) = t \cdot a_0 %2B a_1 x %2B a_2 x^2 %2B \ldots %2B a_{14} x^{14}.">

In this context note that actual Littlewood polynomials is the frame associated to <img src="https://render.githubusercontent.com/render/math?math=t = 1"> and no 
other. Also note that we can't go to infinity, however visual convergence was observed already for t less than 500. 

<img src="https://github.com/felipebottega/PolynomialCodes/blob/master/Littlewood%20Polynomials/im_readme.png"><br><br>

For more about Littlewood polynomials check these links:

https://math.ucr.edu/home/baez/roots/

https://en.wikipedia.org/wiki/Littlewood_polynomial
