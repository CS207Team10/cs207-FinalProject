# Milestone 2

[Milestone 2 Rubric](https://github.com/IACS-CS-207/cs207-F17/blob/master/project/milestone2/milestone2_rubric.md)


## Proposed Feature
You should motivate and describe the feature
Explain how the feature will fit into your code base (and package)
Discuss the modules that you will write to realize your feature
Map out the methods you plan on implementing
Overview how you envision the user to use your new feature
Discuss any external dependencies that your feature will require


Based on Le Chatelier's principle, changing the concentration of a chemical will have an effect to the reaction equilibrium. Thus, the reaction rate, extent, and yield of products will be altered corresponding to the impact on the system. 

For irreversible reactions, users may want to know the progress of the reaction by checking the concentration for each reactant and product. For reversible reactions, users may also want to check the concentration for each specie, in order to determine whether the reaction has reached the equilibrium. 

Our team plan to add a new feature to compute the concentration for each specie at a specific time which will be provided by users. Since the concentration change of specie $i$ is determined by the ODE, 

<img src="http://www.sciweavers.org/tex2img.php?eq=%24%24%5Cfrac%7B%5Cmathrm%7Bd%7Dx_%7Bi%7D%7D%7B%5Cmathrm%7Bd%7Dt%7D%20%3D%20f_%7Bi%7D%5Cleft%28%5Cmathrm%7Bx%7D%2C%20T%5Cright%29%2C%20%5Cqquad%20i%20%3D%201%2C%20%5Cldots%2C%20N.%24%24&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t} = f_{i}\left(\mathrm{x}, T\right), \qquad i = 1, \ldots, N.$$" width="244" height="28" />

and we already had the feature to compute the reaction rates, we can easily get the concentration for each specie at a given time by intergrating the reaction rates by time.
