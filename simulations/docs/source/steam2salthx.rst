.. steam2salthx:

Indirect Configuration Steam-to-Salt Heat Exchanger
######################################################

.. image:: _static/steam2salt_hx.png
   :target: _static/steam2salt_hx.png

Heat exchanger diagram, based on White et al. (2021). Heat lost to surroundings assumed to be 0.

Hot Stream - Steam
-------------------------------------
Hot stream in this case is the fraction of steam coming from the LFR loop.
The heat exchanged between the two streams can be written in terms of the hot stream properties as: 

.. math:: 

   \dot{Q}_{HX} = \dot{m}^H ( h^H_{in} - h^H_{out} ).

For the steam, we have loop-up tables for enthalpy as a function of temperature and pressure. 
We assume that the steam inlet temperature, mass flow and pressure are known:

.. math:: 
   \begin{align}
        \dot{m}^H &= \text{calculated from } \texttt{nuclear_plant} \\
   	T^H_{in} &= 570 {}^\circ \text{C} \\
   	P^H_{in} &= 33 \text{MPa}
   \end{align}

Outlet temperature conditions for the steam are not known ahead of this calculation.

Cold Stream - Molten Salt
-------------------------------------

Cold stream in this case involves the molten salt taken from cold tank and sent to hot tank of thermal energy storage (TES) system. 
The heat exchanged between the two streams can be written in terms of the cold stream properties as:

.. math:: 

   \dot{Q}_{HX} = \dot{m}^C ( h^C_{out} - h^C_{in} ).
   
We want to calculate both the mass flow of the salt. Based on the discussion in `model2mixing`_ we can make an assumption that the change in enthalpy of the molten salt is proportional to the change in temperature

.. math::

   \Delta h^C = c^{\text{salt}}_P ( T^C_{out} - T^C_{in})
   
so we can rewrite the heat exchanged as 

.. math:: 

   \dot{Q}_{HX} = \dot{m}^C c^{\text{salt}}_P (  T^C_{out} - T^C_{in} ).

The known values here are:

.. math:: 
   \begin{align}
   	T^C_{in} &= \text{calculated from } \texttt{tes_two_tank} \\
   	T^C_{out} &= 560 {}^\circ \text{C} \\
   	c^{\text{salt}} &= \text{calculated from } T^C_{in}
   \end{align}
   
Heat Exchanger Efficiency
-------------------------------------

The efficiency of the heat exchanger is written as

.. math::

   \epsilon = \frac{\dot{Q}_{HX}}{\dot{Q}_{max}}

where the maximum heat exchanged is defined under perfect exchange conditions. 
In this case, from the POV of the hot stream we assume that the outlet hot stream temperature is the same as the inlet cold stream temperature:

.. math::

   \dot{Q}_{max} = \texttt{min}(\dot{m}^H, \dot{m}^C) (h^H_{in} - h^C_{in}  ).
   
We then calculate the actual heat exchanged as

.. math::

   \dot{Q}_{HX} = \epsilon \texttt{min}(\dot{m}^H, \dot{m}^C) (h^H_{in} - h^C_{in}  )



Cold Stream Mass Flow
-------------------------------------

The cold stream mass flow can be calculated by combining the two final equations for heat exchanged:

.. math::

   \dot{m}^C = \frac{\epsilon \texttt{min}(\dot{m}^H, \dot{m}^C) (h^H_{in} - h^C_{in}  )}{c^{\text{salt}}_P (  T^C_{out} - T^C_{in} )}


