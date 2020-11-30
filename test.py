import main

# Testing soil weight est
result1 = main.soil_weight_est(1.0, 'kPa', 'psf')
result2 = main.soil_weight_est(.01, 'MPa', 'psf')
result3 = main.soil_weight_est(1.0, 'psi', 'psf')

# Testing Corrected Cone Resistance
result1 = main.corrected_cone_resistance(3000,0)
result2 = main.corrected_cone_resistance(250,10)

# Testing Effective Overburden Pressure
result1 = main.effective_overburden_pressure(117.2,11,above_pressure=1445,depth_below_water_table=6,units='english')
result1 = main.convert(result1,'psf','psi')
here = ""

