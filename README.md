# pyGIM_3DHF

2D and 3D geothermal heat flow, using pyGIMLi finite element solver.

Most geothermal heat flow (GHF) maps for Antarctica are derived in 1D calculations or point-wise manner, use a constant crustal heat production value and disregard the regionally very strong sediment layers. This leads to many errors stemming from lateral effects like heat refraction between geological units and misrepresents the effect of heat production distribution on GHF. This code allows for regional 3D GHF calculations, as well as 2D calculations along profiles.

Input data consists of geophysical boundary layers (topography, Moho-depth, LAB-depth, optional sediment layer-depth). This needs to be in same coordinate system, but does not require to be on same resolution. During domain selection resolution can be chosen accordingly and data will be interpolated on fitting grid structure. Mesh creation uses the pyGIMLi mesh creator, based on some other triangle mesh creator (update here). Alternatively, mesh can be created with external program and given as input for solver.

Boundary conditions of 0° C and 1315° C are applied to the surface and LAB, respectively. Then pyGIMLis finite element sovler is used to derive a temperature field for the model domain. The temperature gradient at surface cells is the corrected by the thermal conductivity of the respective geologic unit to derive surface geothermal heat flow values.

