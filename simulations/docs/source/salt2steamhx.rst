.. _salt2steamhx:

Indirect Salt-to-Steam Heat Exchanger
######################################################

.. image:: _static/salt2steam_hx.png
   :target: _static/salt2steam_hx.png
   :width: 300px
   :align: center

Heat exchanger diagram, based on `White et al. (2021) <https://doi.org/10.3390/su132212428>`_. Heat lost to surroundings assumed to be 0.

Hot Stream - **Molten Salt**
-------------------------------------

Hot stream in this case is the molten salt taken from hot tank and sent to cold tank of thermal energy storage (TES) system. 
The heat exchanged between the two streams can be written in terms of the hot stream properties as: 

.. math:: 
   :name: sEq. 1
   
   \color{red}{\dot{Q}_{HX}} \color{black}{ = \dot{m}^H c_p^H ( T^H_{in} - T^H_{out} ). }

Unknown properties are highlighted in red. The superscript :math:`{}^H` represents the hot stream into the HX which in this case is the **salt**. We assume that the salt inlet temperature, mass flow and outlet temperature are known (calculated within SSC):

.. math:: 
   \begin{align}
        &\dot{m}^H \\
        &c_p^H \\
   	&T^H_{in} \\
   	&T^H_{out}
   \end{align}

From `sEq. 1`_ we can find the total heat exchanged.


Cold Stream - **Steam**
-------------------------------------

Cold stream in this case involves steam through the TES dispatch loop that will recombine downstream with the steam from the LFR loop. 
The heat exchanged between the two streams can be written in terms of the cold stream properties as:

.. math:: 
   :name: sEq. 2
   
   \color{red}{\dot{Q}_{HX}} \color{black}{ = \dot{m}^C ( } \color{red}{h^C_{out}} \color{black}{ - h^C_{in}  ).}
   
The superscript :math:`{}^C` represents the cold stream into the HX which in this case is the **steam**. We want to calculate the mass flow of the steam and outlet conditions. 

The known values here are:

.. math:: 
   \begin{align}
   	T^C_{in} &= 340 {}^\circ \text{C}\\
   	P^C_{in} &= 33 \text{MPa} 
   \end{align}

The enthalpies in `sEq. 2`_ are therefore:

.. math:: 
   \begin{align}
        h^C_{in} &= f(T^C_{in}, P^C_{in}) \\
   	\color{red}{h^C_{out}} &= f(\color{red}{T^C_{out}} \color{black}{ , P^C_{in}) }
   \end{align}
   
Heat Exchanger Effectiveness
-------------------------------------

The effectiveness of the heat exchanger is written as

.. math::
   :name: sEq. 3
   
   \color{black}{\epsilon = } \frac{\color{red}{\dot{Q}_{HX}}}{\color{red}{\dot{Q}_{max}}}

where the maximum heat exchanged is defined under perfect heat exchange conditions. 
In this case, we use the minimum capacitance rate of the two fluids (assuming cold fluid has lower minimum capacitance rate). We also use the maximum difference in enthalpies:

.. math::
   :name: sEq. 4
   
   \color{red}{\dot{Q}_{max}} \color{black}{ = \dot{m}^C ( h^C_{out,ideal} - h^C_{in}  ).}

where

.. math::

   h^C_{out,ideal} = f(T^H_{in}, P^C_{in})

which implies that under these ideal conditions, the steam will exit the HX with the maximum temperature change. We then calculate the actual heat exchanged by combining `sEq. 3`_ and `sEq. 4`_ as

.. math::
   :name: sEq. 5
   
   \color{red}{\dot{Q}_{HX}} \color{black}{ = \epsilon \dot{m}^C (h^C_{out,ideal} - h^C_{in}  ). }



Solving for Unknown Variables
-------------------------------------

The cold stream mass flow can be calculated by combining the `sEq. 1`_ and `sEq. 5`_ as:

.. math::

   \color{red}{\dot{m}^C} \color{black}{ = \frac{ \dot{m}^H c^{\text{salt}}_P (  T^C_{in} - T^C_{out} )}{\epsilon \dot{m}^H (h^H_{out,ideal} -h^H_{in} )} }


The outlet temperature of the steam within the HX can be backed out by solving for :math:`\color{red}{\dot{Q}_{HX}}` in `sEq. 1`_ and substituting that into `sEq. 2`_ to find:

.. math:: 
   
   \color{red}{h^C_{out}} \color{black}{ = h^C_{in} + \frac{\dot{Q}_{HX}}{\dot{m}^C } }

from which we can back out the temperature by through the steam tables

.. math::

   \color{red}{T^C_{out}} \color{black}{ = f^{-1}(h^C_{out}, P^C_{in} ) } 

