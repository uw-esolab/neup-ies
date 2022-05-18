.. designPointMassFlow:

Design Point Mass Flow Calculations
#################################################


Calculations in the Direct Charging Case
--------------------------------------------

In the direct TES charging design (e.g., :math:`\texttt{cmod_tcsmolten_salt}, \texttt{cmod_nuclear_tes}, \texttt{cmod_nuclear_mspt_tes}`) we create a power cycle object. 

The following happens in the :math:`\texttt{init()}` method of the :math:`\texttt{csp_solver_pc_Rankine_indirect_224.cpp}` class.

Using inputs from a JSON script, we can specify the reference power output :math:`\texttt{P_ref}` and a design point efficiency :math:`\texttt{eta}`. From these two parameters, we calculate a design point thermal power input :math:`\texttt{m_q_dot_design}`. 

.. math::

   \texttt{m_q_dot_design} = \frac{\texttt{P_ref}}{\texttt{eta}}

From there, the design point mass flow :math:`\texttt{m_q_dot_design}` is calculated from

.. math::

   \texttt{m_m_dot_design} = \frac{\texttt{m_q_dot_design}}{\texttt{m_cp_htf_design} * (\texttt{m_T_htf_hot_ref} - \texttt{m_T_htf_cold_ref})}

and converted to the appropriate units. The specific heat here :math:`\texttt{m_cp_htf_design}` is of the heat transfer fluid which, in the nominal case, is #17 in the enumeration which corresponds to molten salt (*Salt_60_NaNO3_40_KNO3*). The molten salt is calculated using empirical formulas from the :math:`\texttt{HTFProperties}` class as a function of temperature:

.. math::

   \texttt{m_cp_htf_design} = f( \frac{1}{2} (\texttt{m_T_htf_hot_ref} + \texttt{m_T_htf_cold_ref} ) ).


This formula holds for the HTF fluid because the salt specific heat does not vary *too* much between those temperature ranges. 

Typical temperatures used here are

.. math::
   \begin{align}
      \texttt{m_T_htf_hot_ref}  &= 560 {}^\circ \text{C} \\
      \texttt{m_T_htf_cold_ref} &= 340 {}^\circ \text{C}
   \end{align}

Off-design behavior for this configuration is determined by looking at the **salt** mass flow through the salt-to-steam HX: values are taken from a user-defined table which should capture efficiency hits and other parameters as a function of hot temperature and mass flow as a fraction of :math:`\texttt{m_m_dot_design}`.


Calculations for the Indirect Charging Case
----------------------------------------------

For the indirect charging case, we care more about the mass flow of the power cycle working fluid rather than the salt. 

So instead of salt, we're looking at steam. The previous equation for mass flow no longer holds because the specific heat of steam varies by a lot between those temperature ranges. Instead we need to use enthalpies.

We can imagine a control volume before and after the heat addition for this configuration (ignoring complexities of stream splitting for LFR and TES loops). The conservation of energy equation looks like:

.. math::

   \dot{q}^{\text{des}} = \dot{m}^{\text{s}} [ \Delta h^{\text{s}}_{H\rightarrow\text{sat}} + \Delta h^{\text{s}}_{sat} + \Delta h^{\text{w}}_{\text{sat}\rightarrow C} ]
 
where the superscripts :math:`{}^{\text{s}}` and :math:`{}^{\text{w}}` are supercritical steam and liquid water, respectively. 

Critical temperature assumed to be :math:`374 {}^\circ` C for pressures above critical pressure 22 MPa. Actual pressure of steam here is 33 MPa. 

In SCC, we can calculate all the terms using the :math:`\texttt{water_properties}` class.
The enthalpy change from hot temperature to cold temperature is:

.. math::

   \begin{align}
      &\texttt{wp}.\texttt{water_TP}( \texttt{m_T_htf_hot_ref} , \texttt{P_boil} )    \\
      &\texttt{h_steam_hot} = \texttt{wp}.\texttt{enth} \\
      &\texttt{wp}.\texttt{water_TP}( \texttt{m_T_htf_cold_ref} , \texttt{P_boil} )    \\
      &\texttt{h_steam_cold} = \texttt{wp}.\texttt{enth} \\
      &\texttt{delta_h} =   \texttt{h_steam_hot} - \texttt{h_steam_cold}.
   \end{align}
   
Because we are operating the steam above the critical pressure of 22 MPa, there is no latent heat of vaporization, so 

.. math::

   \Delta h^{\text{s}}_{sat} = 0

and the energy conservation equation becomes:

.. math::

   \dot{q}^{\text{des}} = \dot{m}^{\text{s}} \Delta h^{\text{s}}_{H\rightarrow C}.
 
Design point mass flow is therefore

.. math::

   \texttt{m_m_dot_design} = \frac{\texttt{m_q_dot_design}}{\texttt{delta_h}}

