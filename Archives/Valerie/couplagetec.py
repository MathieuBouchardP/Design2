from scipy.optimize import curve_fit

tension = [-0.3, -0.3, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
puissance = [0.945, 0.936, 0.63, 0, 0.67, 1.78, 2.2275, 2.82]

def couplage_function(x, a, b, c):
  return a*x**2 + b*x + c

# curve_fit
popt, pcov = curve_fit(f=couplage_function, xdata=tension, ydata=puissance)
print(popt)